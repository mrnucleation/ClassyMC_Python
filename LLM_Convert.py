import os
import subprocess
import torch

from typing import List, Optional

import fire
import sys


#sys.path.append('../Llama/')

from llama import Dialog, Llama


# ------------------------------------------------------------------------------
def run_llama(
    prompt: str,
    ckpt_dir: str,
    tokenizer_path: str,
    temperature: float = 0.6,
    top_p: float = 0.9,
    max_seq_len: int = 3012,
    max_batch_size: int = 4,
    max_gen_len: Optional[int] = None,
):
    """
    Examples to run with the models finetuned for chat. Prompts correspond of chat
    turns between the user and assistant with the final one always being the user.

    An optional system prompt at the beginning to control how the model should respond
    is also supported.

    The context window of llama3 models is 8192 tokens, so `max_seq_len` needs to be <= 8192.

    `max_gen_len` is optional because finetuned models are able to stop generations naturally.
    """
    generator = Llama.build(
        ckpt_dir=ckpt_dir,
        tokenizer_path=tokenizer_path,
        max_seq_len=max_seq_len,
        max_batch_size=max_batch_size,
    )

    dialogs: List[Dialog] = [
        [{"role": "user", "content": prompt}]]
    
    results = generator.chat_completion(
        dialogs,
        max_gen_len=max_gen_len,
        temperature=temperature,
        top_p=top_p,
    )

    
    for dialog, result in zip(dialogs, results):
        for msg in dialog:
            print(f"{msg['role'].capitalize()}: {msg['content']}\n")
        print(
            f"> {result['generation']['role'].capitalize()}: {result['generation']['content']}"
        )
        print("\n==================================\n")
    return results[0]
# ------------------------------------------------------------------------------

# Check if CUDA is available
if torch.cuda.is_available():
    device = torch.device("cuda")
    print("Using CUDA")
else:
    raise Exception("CUDA is not available. Please make sure CUDA is installed and configured.")


# Specify the directory containing the Fortran files

workdir = os.getcwd()
fortran_files_dir = workdir+'/src/'

from glob import glob
filenames = glob(fortran_files_dir + '*.f90')
filenames = sorted(filenames)

filenames = filenames[:1]

# Change to the directory containing the Fortran files
os.chdir(fortran_files_dir)

prompt = "Take this fortran 2003 style code and convert it into Python 3.  For any arrays make sure to use numpy. Any classes should be converted faithfully."
prompt += " Please also make sure to convert comments over to Python 3 style comments. "

# Loop through all files in the directory
for filename in filenames:
    outfile = filename.replace('.f90', '.py')
    # Check if the file has a .f90 extension
    # Load the contents of the f90 file
    with open(filename, 'r') as f:
        f90_code = f.read()
    
    # Append the f90 code to the prompt string
    prompt += f"\n\n{f90_code}"
    
    print(prompt)
    
    result = fire.Fire(run_llama, args=[prompt])
    