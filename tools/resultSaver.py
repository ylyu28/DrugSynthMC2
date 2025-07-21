from pathlib import Path
import json

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



def load_update_mast(filename, ngram_list, kd):
    try:
        with open(filename, 'r') as f:
            mast_dict = json.load(f)

        if ngram_list is not None and kd is not None:
            for ngram in ngram_list: 
                if ngram in mast_dict:
                    mast_dict[ngram].append(kd)
                else:
                    mast_dict[ngram] = [kd]

        with open(filename, 'w') as f:
            json.dump(mast_dict, f, indent=4)
    
    except FileNotFoundError:
        mast_dict = {}
        if ngram_list is not None and kd is not None:
            for ngram in ngram_list:
                mast_dict[ngram] = [kd]
        with open(filename, 'w') as f:
            json.dump(mast_dict, f, indent=4)
    
    except Exception as e:
        print(f"Error occurred:{e}")



def calc_mast_ngram(mast_ngram_len, smiles):

    s = ''
    for i in range(len(smiles)):
        c = smiles[i]
        # substitute ring opening with X
        if c.isdigit() and smiles[:i+1].count(c) == 1:
            s += 'X'
        else:
            s += c
    print(s)
    ngram_s = list(set(s[i:mast_ngram_len+i] for i in range(len(s) - mast_ngram_len)))
    
    # ignore the mast ngrams that contain ring closure since this move is set fixed probabilities
    filtered_mast_ngrams = [ngram for ngram in ngram_s if not any(char.isdigit() for char in ngram)]

    return filtered_mast_ngrams