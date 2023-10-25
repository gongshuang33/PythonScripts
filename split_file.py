import argparse


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-f', type=str, required=True, help='Input one file.(STR, required)')
    parser.add_argument('--part', '-p', type=int, required=True, help='Slice the specified number of files.(INT, required)')
    args = parser.parse_args()
    return args

def main():
    args()

if __name__ == '__main__':
    main()