import psutil
import time

# Configuration parameters
THRESHOLD_MB = 62700  # Memory usage threshold in MB for 4x64x64 Offline jobs on Milan
CHECK_INTERVAL = 10  # Check interval in seconds

def get_processes_by_name(process_name):
    """Finds all processes by name."""
    matching_processes = []
    for process in psutil.process_iter(['name']):
        if process.info['name'] == process_name:
            matching_processes.append(process)
    return matching_processes

def check_memory_usage_and_kill(process_name, threshold_mb):
    """Checks the memory usage of all processes with a given name and kills them if they exceed the threshold."""
    processes = get_processes_by_name(process_name)

    if not processes:
        print(f"No processes named '{process_name}' found.")
        return

    for process in processes:
        try:
            # Get the memory usage in MB
            memory_usage_mb = process.memory_info().rss / (1024 * 1024)
            print(f"Memory usage of '{process_name}' (PID: {process.pid}): {memory_usage_mb:.2f} MB vs threshold {THRESHOLD_MB} MB")

            # Check if the memory usage exceeds the threshold
            if memory_usage_mb > threshold_mb:
                print(f"Memory usage exceeded threshold ({threshold_mb} MB). Killing process {process.pid}.")
                process.terminate()  # Send a termination signal
                try:
                    process.wait(timeout=60)  # Wait for the process to terminate, in seconds
                except psutil.TimeoutExpired as e:
                    print(f"Error wating for process {process.pid} to terminate: {e}")
                print(f"Process {process.pid} terminated.")
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess) as e:
            print(f"Error accessing process {process.pid}: {e}")

if __name__ == "__main__":
    while True:
        check_memory_usage_and_kill('cmsRun', THRESHOLD_MB)
        time.sleep(CHECK_INTERVAL)

