#!/usr/bin/env python3

import wehrli_gauss

def main():
    wehrli_gauss.filter_file("data/wehrli_original.csv", "data/wehrli_filtered.csv")

if __name__ == "__main__":
    main()