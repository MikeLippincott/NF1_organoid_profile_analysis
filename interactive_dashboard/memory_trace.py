"""Process memory tracing for the dashboard.

Reads the resident set size (RSS) of the running Python process. On Linux that
comes from ``/proc/self/status``; elsewhere it falls back to ``psutil`` if it is
installed. The peak comes from ``resource`` on every platform.

- ``trace_memory(label)`` logs the RSS change across a block to stderr.
- ``sidebar_memory()`` shows the current and peak RSS in the sidebar.
"""

import contextlib
import resource
import sys

import streamlit as st


def current_rss_mb() -> float:
    """Current resident memory in MB, or NaN if it can't be read."""
    try:
        with open("/proc/self/status") as status:
            for line in status:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024
    except OSError:
        pass
    try:
        import psutil

        return psutil.Process().memory_info().rss / 1024**2
    except ImportError:
        return float("nan")


def peak_rss_mb() -> float:
    """Highest resident memory this process has reached, in MB."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is KB on Linux and bytes on macOS
    return peak / 1024**2 if sys.platform == "darwin" else peak / 1024


@contextlib.contextmanager
def trace_memory(label: str):
    """Print the RSS change across the block to stderr (also on error)."""
    before = current_rss_mb()
    try:
        yield
    finally:
        after = current_rss_mb()
        print(
            f"[memory] {label}: {before:,.0f} -> {after:,.0f} MB "
            f"(delta {after - before:+,.0f}, peak {peak_rss_mb():,.0f})",
            file=sys.stderr,
            flush=True,
        )


def sidebar_memory() -> None:
    """Live memory readout, refreshed on every rerun."""
    st.sidebar.subheader("Memory")
    current = current_rss_mb()
    st.sidebar.metric(
        "Process RSS",
        f"{current:,.0f} MB" if current == current else "unavailable",
        help="Resident memory of this app process, read on every rerun.",
    )
    st.sidebar.caption(f"Peak since start: {peak_rss_mb():,.0f} MB")
