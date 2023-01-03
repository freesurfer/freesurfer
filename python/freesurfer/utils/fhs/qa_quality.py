#!/usr/bin/env python3
"""
USAGE: ./tp_quality </path/to/mancheck.json> [--require_second_qa (optional) will check for 'second_quality_check']

Checks to see if the file exists, if so, prints the value of 'overall_score' to terminal
If file does not exist or 'overall_score' is not in file, prints -1
If file cannot be accessed, prints -2
If 'second_quality_check' is missing or false, prints -3

J. Nolan 09/21/21
Edited 10/26/21 J. Nolan
    - Added check for 'second_quality_check'
    - If 'second_quality_check' not present, returns -3
Edited 11/4/21 J. Nolan
    - Added flag to make "second_quality_check" an optional check
    - Addressed issue validating "second_quality_check" value
"""
import os, json, argparse

def parse_args():
    parser = argparse.ArgumentParser(usage='./tp_quality </path/to/mancheck.json>\n\
        Prints the value of "overall_score" from mancheck.json to terminal')
    parser.add_argument('path',type=str, help='Path to mancheck.json')
    parser.add_argument('--require_second_qa',help='(optional) checks for "second_quality_check" in the mancheck.json', action='store_true',default=False)

    return parser.parse_args()

def get_overall_score(path_to_mancheck,second_qa=False):
    if not (os.path.isfile(path_to_mancheck)):
        return -1
    try:
        f = open(path_to_mancheck)
        data = json.load(f)
        f.close()

        if('overall_score' not in data.keys()):
            return -1

        quality = data['overall_score']

        if(second_qa):
            if ('second_quality_check' in data.keys()):
                if (data['second_quality_check'] == True):
                    return quality
            elif ('qa_level' in data.keys()):
                if (data['qa_level'] == 4):
                    return quality
            return -3

        return quality

    except:
        return -2

def main():
    args = parse_args()
    print(get_overall_score(args.path,args.require_second_qa))

if __name__ == '__main__':
    main()