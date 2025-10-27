import sys

def print_progress(t, start, duration, BAR_LENGTH = 30):
    """Progress bar for solvers."""

    last_update = -1
    progress = int(100 * (t - start) / duration)
    if progress > last_update:
        filled = int(BAR_LENGTH * progress / 100)
        bar = "█" * filled + "-" * (BAR_LENGTH - filled)
        sys.stdout.write(f"\rProgress: |{bar}| {progress:3d}%")
        sys.stdout.flush()
        last_update = progress