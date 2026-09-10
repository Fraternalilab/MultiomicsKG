import csv
import time
import datetime
import gzip
import os

NUM_ITEMS_PER_PROGRESS_REPORT = 10000

class ProgressFile:

    def __init__(self, file, compressed=False):

        self.compressed = compressed
        self.file = file
        if self.compressed:
            self.raw = self.file.fileobj.raw
        self.size = os.fstat(file.fileno()).st_size

    def __iter__(self):
        return self

    def tell(self):
        if self.compressed:
            return self.raw.tell()
        else:
            return self.file.tell()

    def __next__(self):
        line = self.file.readline()

        if not line:
            raise StopIteration
        return line.decode("utf-8")


class ICSV:
    """Represents a generic CSV file interface"""

    def __init__(
        self,
        filename: str,
        file_type: str,
        max_items=-1,
        header_filter_prefix="",
        delimiter=",",
        quotechar='"',
        skip_items=0,
        compressed=False,
        skip_rows=0,
        fieldnames=None,
    ):
        self.file_type = file_type
        self.skip_items = skip_items
        self.max_items = max_items
        self.num_imported = 0
        self.filename = filename
        self.skip_rows = skip_rows
        if compressed:
            self.filestream = gzip.open(filename=filename, mode="rb")
            self.source_file = ProgressFile(self.filestream, compressed=True)
        else:
            self.filestream = open(file=filename, mode="rb")
            self.source_file = ProgressFile(self.filestream, compressed=False)
        if self.skip_rows > 0:
            i = 0
            while i < self.skip_rows:
                next(self.source_file, None)
                i += 1
            self.reader = csv.DictReader(
                self.source_file,
                delimiter=delimiter,
                quotechar=quotechar,
                fieldnames=fieldnames,
            )
        elif len(header_filter_prefix) > 0:
            self.reader = csv.DictReader(
                filter(
                    lambda row: not row.startswith(header_filter_prefix),
                    self.source_file,
                ),
                delimiter=delimiter,
                quotechar=quotechar,
                fieldnames=fieldnames,
            )
        else:
            self.reader = csv.DictReader(
                self.source_file,
                delimiter=delimiter,
                quotechar=quotechar,
                fieldnames=fieldnames,
            )
        # self.filestream.seek(0, os.SEEK_END)
        # self.filestream_size = self.filestream.tell()
        # self.filestream.seek(0, os.SEEK_SET)
        self.items_read = 0
        self.report_line = (
            self.file_type
            + " items imported: {num_items:7d}   |   "
            + "Progress: {progress:8.2f} %  |   "
            + "Time elapsed: {time_elapsed}   |   "
            + "ETA: {eta}               "
        )

    def get_dict_items(self, user_key: str) -> dict():
        # need completion with counters and skip items
        data = dict()
        while True:
            item = next(self.reader, None)
            if item is None:
                return data

            self.items_read += 1
            key = item[user_key]

            if key in data:
                data[key].append(item)
            else:
                data[key] = list()
                data[key].append(item)

    def importfile(self):
        start_counter = time.perf_counter()
        stop_counter = start_counter
        stop_position = start_position = 0
        start = time.time()

        # skip items
        while True:
            while self.items_read < self.skip_items:
                if self.items_read == self.max_items:
                    return
                item = next(self.reader, None)
                self.items_read += 1
                continue

            if self.items_read == self.max_items:
                return

            item = next(self.reader, None)

            if item is None:
                time_elapsed = str(
                    datetime.timedelta(seconds=round(time.time() - start))
                )
                print(
                    self.report_line.format(
                        num_items=self.num_imported,
                        progress=100.0,
                        time_elapsed=time_elapsed,
                        eta=0,
                    ),
                    end="\n",
                )
                return

            self.items_read += 1

            self.process_item(item)

            self.num_imported += 1
            if self.items_read % NUM_ITEMS_PER_PROGRESS_REPORT == 0:
                progress = (self.source_file.tell() / self.source_file.size) * 100.0
                stop_counter = time.perf_counter()
                stop_position = self.source_file.tell()
                speed = (stop_position - start_position) / (
                    stop_counter - start_counter
                )
                eta = str(
                    datetime.timedelta(
                        seconds=round((self.source_file.size - stop_position) / speed)
                    )
                )
                time_elapsed = str(
                    datetime.timedelta(seconds=round(time.time() - start))
                )
                print(
                    self.report_line.format(
                        num_items=self.num_imported,
                        progress=progress,
                        time_elapsed=time_elapsed,
                        eta=eta,
                    ),
                    end="\r",
                )
                start_counter = time.perf_counter()
                start_position = stop_position

    def close(self):
        self.filestream.close()

    def get_next_item(self):
        if self.items_read == self.max_items:
            return None

        item = next(self.reader, None)

        if item is None:
            print(
                "{file_type} items read: {num_items:7d}, progress: {progress:8.2f} %".format(
                    file_type=self.file_type,
                    num_items=self.items_read,
                    progress=100.0,
                ),
                end="\n",
            )
            return None

        self.items_read += 1
        if self.items_read % NUM_ITEMS_PER_PROGRESS_REPORT == 0:
            progress = (self.source_file.tell() / self.source_file.size) * 100.0
            print(
                "{file_type} items read: {num_items:7d}, progress: {progress:8.2f} %".format(
                    file_type=self.file_type,
                    num_items=self.items_read,
                    progress=progress,
                ),
                end="\r",
            )

        return item
