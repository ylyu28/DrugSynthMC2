from pathlib import Path

def writeline(line:str, name:str):
    if name:
        file_path = Path('result')/f"{name}.csv"
    try:
        file_path.parent.mkdir(parents=True, exist_ok=True)
        if file_path.exists():
            with file_path.open("a") as file:
                file.write(line)
        else:
            with file_path.open("w") as file:
                file.write(line)
    except (OSError, IOError) as e:
        print(f"Error writing to file: {e}")
