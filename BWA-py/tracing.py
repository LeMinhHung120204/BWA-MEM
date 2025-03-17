from dataclasses import dataclass, asdict
from typing import Literal, List
import time
import os
import ctypes
import threading
import json
import shutil
from pathlib import Path
import uuid

libc = ctypes.cdll.LoadLibrary("libc.so.6")


def _get_tid() -> int:
    return threading.get_ident()

    SYS_gettid = 186
    tid = libc.syscall(SYS_gettid)
    return tid


@dataclass
class TraceEvent:
    # https://docs.google.com/document/d/1CvAClvFfyA5R-PhYUmn5OOQtYMH4h6I0nSsKchNAySU/preview?tab=t.0
    # all time is in microsecond here
    name: str
    cat: str
    ph: Literal["B", "E", "X"]
    pid: int
    tid: int
    ts: float
    dur: float


class TraceContext:
    def __init__(self, tracer: "Tracer", name: str):
        self.tracer = tracer
        self._name = name
        self._ts = None

    def __enter__(self) -> "TraceContext":
        self._ts = self.tracer.get_timestamp()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            return False

        self.tracer.record(self._name, self._ts)
        return False


class Tracer:
    def __init__(self, trace_dir: str):
        # self._pid = os.getpid()
        # self._tid = _get_tid()

        # perfectto UI group threads base on pid, to make them appear nicely, we
        # set them to have the same pid while having different thread id.
        self._pid = 1
        self._tid = os.getpid()
        self._uid = uuid.uuid4().hex

        self._trace_file = (
            Path(trace_dir) / f"{self._pid}_{self._tid}_{self._uid}.trace"
        )
        self._trace_file.touch()

    def get_timestamp(self) -> int:
        return time.perf_counter_ns()

    def record(self, event_name: str, begin_timestamp_ns: int) -> TraceEvent:
        dur_ns = time.perf_counter_ns() - begin_timestamp_ns
        event = TraceEvent(
            name=event_name,
            cat="custom",
            ph="X",
            pid=self._pid,
            tid=self._tid,
            ts=begin_timestamp_ns / 1e3,
            dur=dur_ns / 1e3,
        )

        with self._trace_file.open("a") as f:
            json.dump(asdict(event), f)
            f.write(",\n")

        return event

    def profile(self, event_name: str) -> TraceContext:
        return TraceContext(self, event_name)


def clean_trace(trace_dir: str):
    shutil.rmtree(trace_dir, ignore_errors=True)
    Path(trace_dir).mkdir()


def make_chrome_trace(trace_dir: str, filename: str = "trace.json") -> Path:
    def is_trace(p: Path):
        return p.suffix == ".trace"

    p = Path(trace_dir)
    files = list(filter(is_trace, p.iterdir()))
    tc = "["
    for file in files:
        with file.open() as f:
            tc += f.read()
    tc += "]"

    trace_file = p / filename
    with trace_file.open("w+") as f:
        f.write(tc)

    return trace_file


if _name_ == "_main_":
    trace_dir = "C:\\Hung\\NCKH\\BWA-py\\trace"

    clean_trace(trace_dir)

    tracer = Tracer(trace_dir)
    with tracer.profile("sleep event"):
        time.sleep(3)

    with tracer.profile("other sleep event"):
        time.sleep(1)

    ts = tracer.get_timestamp()
    time.sleep(2)
    tracer.record("manually record event", ts)

    trace_file = make_chrome_trace(trace_dir)
    print("collected chrome trace file at:", trace_file)