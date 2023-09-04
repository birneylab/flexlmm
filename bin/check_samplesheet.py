#!/usr/bin/env python


"""
Provide a command line tool to validate and transform tabular samplesheets.
"""


import argparse
import csv
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a
            previously validated and transformed row.
            The order of rows is maintained.
    """

    VALID_CRAM_FORMATS = [".cram"]
    VALID_CRAI_FORMATS = [".crai"]
    VALID_HIGH_COV_VALUES = ["true", "false"]

    def __init__(
        self,
        sample_col="sample",
        first_col="cram",
        second_col="crai",
        third_col="high_cov",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample
                name (default "sample").
            first_col (str): The name of the column that contains the cram file
                path (default "cram").
            second_col (str): The name of the column that contains the crai
                file (default "crai").
            third_col (str): The name of the new column that will be inserted
                if not present and records whether a sample is high coverage
                and is to be used in validation.

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._third_col = third_col
        self._seen = {
            self._sample_col: set(),
            self._first_col: set(),
            self._second_col: set(),
        }
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the high_cov column
        if absent.

        Args:
            row (dict): A mapping from column headers (keys) to elements of
                that row (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_third(row)
        self._seen[self._sample_col].add(row[self._sample_col])
        self._seen[self._first_col].add(row[self._first_col])
        self._seen[self._second_col].add(row[self._second_col])
        self.modified.append(row)

    def _validate_sample(self, row):
        """
        Assert that the sample name exists and convert spaces to underscores.
        """

        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the cram entry is non-empty and has the right format."""

        if len(row[self._first_col]) <= 0:
            raise AssertionError("The cram file is required.")
        self._validate_cram_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the crai entry is non-empty and has the right format"""

        if len(row[self._first_col]) <= 0:
            raise AssertionError("The crai file is required.")
        self._validate_crai_format(row[self._second_col])

    def _validate_third(self, row):
        """
        Assert that the high_cov entry contains only boolean values if present,
        and initialise it to False if absent
        """

        if self._third_col not in row:
            row[self._third_col] = "false"

        self._validate_high_cov_format(row[self._third_col])

    def _validate_cram_format(self, filename):
        """
        Assert that a given filename has one of the expected cram extensions.
        """

        if not any(
            filename.endswith(extension)
            for extension in self.VALID_CRAM_FORMATS
        ):
            raise AssertionError(
                f"The cram file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_CRAM_FORMATS)}"
            )

    def _validate_crai_format(self, filename):
        """
        Assert that a given filename has one of the expected crai extensions.
        """

        if not any(
            filename.endswith(extension)
            for extension in self.VALID_CRAI_FORMATS
        ):
            raise AssertionError(
                f"The crai file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_CRAI_FORMATS)}"
            )

    def _validate_high_cov_format(self, value):
        """
        Assert that a given value is "true" or "false".
        """

        if not any(
            value == extension for extension in self.VALID_HIGH_COV_VALUES
        ):
            raise AssertionError(
                f"The high_cov field has an unrecognized value: {value}\n"
                f"It should be one of: {', '.join(self.VALID_HIGH_COV_VALUES)}"
            )

    def validate_unique_samples(self):
        """
        Assert that the sample name, cram file, and crai file are all unique.
        """

        if len(self._seen[self._sample_col]) != len(self.modified):
            raise AssertionError("Sample names must be unique.")

        if len(self._seen[self._first_col]) != len(self.modified):
            raise AssertionError("Cram filenames must be unique.")

        if len(self._seen[self._second_col]) != len(self.modified):
            raise AssertionError("Crai filenames must be unique.")


def read_head(handle, num_lines=10):
    """
    Read the specified number of lines from the current position in the file.
    """
    lines = []

    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)

    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read
        position is expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)

    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core
    pipelines.

    Validate the general shape of the table, expected columns, and each row. If
    a column called high_cov is not present, it is created and assigned a value
    of False.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can
            be either CSV, TSV, or any other format automatically recognized by
            ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed
            samplesheet should be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following
        structure

            sample,cram,crai
            SAMPLE,SAMPLE.cram,SAMPLE.cram.crai
            SAMPLE,SAMPLE.cram,SAMPLE.cram.crai
    """
    required_columns = {"sample", "cram", "crai"}
    optional_columns = {"high_cov"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on
    # `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.

        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(
                "The sample sheet **must** contain these column headers:"
                + f" {req_cols}."
            )
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()

        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} on line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(required_columns | optional_columns)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on
    # `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()

        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog=(
            "Example: python check_samplesheet.py samplesheet.csv"
            + " samplesheet.valid.csv"
        ),
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )

    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level, format="[%(levelname)s] %(message)s"
    )

    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
