import sys
import time
import logging
from watchdog.observers import Observer
from watchdog.events import LoggingEventHandler, FileSystemEventHandler
from collections import defaultdict


class OnMyWatch:

    def __init__(self, watchDirectory):
        self.observer = Observer()
        self.watchDirectory = watchDirectory

    def run(self):
        event_handler = Handler()
        self.observer.schedule(
            event_handler, self.watchDirectory, recursive=True)
        self.observer.start()
        try:
            while True:
                time.sleep(5)
        except:
            self.observer.stop()
            logging.info("Observer Stopped")

        self.observer.join()


class Handler(FileSystemEventHandler):

    counter = defaultdict(int)

    @classmethod
    def on_any_event(cls, event):
        if event.is_directory:
            return None
        else:
            cls.counter[event.event_type] += 1
            logging.info(", ".join(f"{num} files has been {oper}" for oper,
                         num in cls.counter.items()))


if __name__ == "__main__":
    if __name__ == '__main__':
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
        path = sys.argv[1] if len(sys.argv) > 1 else '.'
        mode = sys.argv[2] if len(sys.argv) > 2 else 'default'
        if mode == 'default':
            event_handler = LoggingEventHandler()
            observer = Observer()
            observer.schedule(event_handler, path, recursive=True)
            observer.start()
            try:
                while True:
                    time.sleep(5)
            except KeyboardInterrupt:
                observer.stop()
            observer.join()
        else:
            watch = OnMyWatch(path)
            watch.run()
