# testing remote exec
import torch
if __name__ == "__main__":
    print(f"is cuda available:{torch.cuda.is_available()}")