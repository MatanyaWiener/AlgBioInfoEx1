import argparse
import numpy as np
from itertools import groupby
from Bio import SeqIO
import csv

from enum import Enum


class Scores(Enum):
    GAP = 0
    A = 1
    C = 2
    G = 3
    T = 4


def print_table(table):
    for row in (table):
        print(row)
    print("\n")

def print_results(seq1, seq2, type, score):
    i = 0 # how many tims we have print 50 letter
    limit = 50
    length = len(seq1)
    while length > i * limit:
        print(seq1[i * limit:min(((i + 1) * limit), length)])
        print(seq2[i * limit:min(((i + 1) * limit), length)], "\n")
        i += 1
    print(type + ":" + str(score))



def fastaread(fasta_name):
    """
    Read a fasta file. For each sequence in the file, yield the header and the actual sequence.
    In Ex1 you may assume the fasta files contain only one sequence.
    You may keep this function, edit it, or delete it and implement your own reader.
    """
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def read_fasta_file(file_path):
    sequences = []
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
    return sequences



def read_tsv_file(file_path):
    data = []
    with open(file_path, 'r', newline='', encoding='utf-8') as tsv_file:
        # Create a CSV reader with tab as the delimiter
        tsv_reader = csv.reader(tsv_file, delimiter='\t')

        # Iterate over rows in the TSV file
        for row in tsv_reader:
            data.append(row)

    return data

def get_score_dictionary(s_mat):
    A = 1
    C = 2
    G = 3
    T = 4
    GAP = 5

    res = {
        ("A", "A"): s_mat[A][A], ("C", "A"): s_mat[A][C], ("G", "A"): s_mat[A][G], ("T", "A"): s_mat[A][T],

        ("A", "C"): s_mat[C][A], ("C", "C"): s_mat[C][C], ("G", "C"): s_mat[C][G], ("T", "C"): s_mat[C][T],

        ("A", "G"): s_mat[G][A], ("C", "G"): s_mat[G][C], ("G", "G"): s_mat[G][G], ("T", "G"): s_mat[G][T],

        ("A", "T"): s_mat[T][A], ("C", "T"): s_mat[T][C], ("G", "T"): s_mat[T][G], ("T", "T"): s_mat[T][T],

        ("A", "-"): s_mat[GAP][A], ("C", "-"): s_mat[GAP][C], ("G", "-"): s_mat[GAP][G], ("T", "-"): s_mat[GAP][T],

        ("-", "A"): s_mat[GAP][A], ("-", "C"): s_mat[GAP][C], ("-", "G"): s_mat[GAP][G], ("-", "T"): s_mat[GAP][T],

    }
    return res


def cell_score__global_matanya(seq1, seq2, score_dictionary, r, c, score_table, path_table):
    """
    We will assume that cell [0][0] is already filled up for us to avoid an "if" statement.
    and that the score matrix is initiated with -inf
    """
    if r == 0 and c == 0:
        score_table[0][0] = 0
        path_table[0][0] = "start"
        return
    # Left check
    letter_col = None
    letter_row = None
    if c > 0:
        letter_col = seq1[c - 1]
        # todo see why in the big function Elisheva used int()
        score_table[r][c] = score_table[r][c-1] + int(score_dictionary.get((letter_col, "-")))
        path_table[r][c] = (letter_col, "-", "left")
    # Up check
    if r > 0:
        letter_row = seq2[r - 1]
        up_score = score_table[r-1][c] + int(score_dictionary.get(("-", letter_row)))
        if score_table[r][c] < up_score:
            score_table[r][c] = up_score
            path_table[r][c] = ("-", letter_row, "up")

    if r > 0 and c > 0:
        up_left_score = score_table[r-1][c-1] + int(score_dictionary.get((letter_col, letter_row)))
        if up_left_score >= score_table[r][c]:
            score_table[r][c] = up_left_score
            path_table[r][c] = (letter_col, letter_row, "up_and_left")


def calculate_score_and_path_matanya(seq1, seq2, score_dictionary, type):
    score_table = [[(float("-inf")) for c in range(len(seq1) + 1)] for r in range(len(seq2) + 1)]
    path_table = [[None for c in range(len(seq1) + 1)] for r in range(len(seq2) + 1)]
    for r in range(len(seq2) + 1):
        for c in range(len(seq1) + 1):
            cell_score__global_matanya(seq1, seq2, score_dictionary, r, c, score_table, path_table)
    return score_table[-1][-1], path_table



def calculate_score_and_path(seq1, seq2, score_dictionary, type):
    result_table = [[(0) for j in range(len(seq1) + 1)] for i in range(len(seq2) + 1)]
    map_matrix = [[("","", "") for j in range(len(seq1) + 1)] for i in range(len(seq2) + 1)]
    max_cell = (0, 0, 0)  # represent the max score in "local" case and the cell's indexes.
    for i in range(len(seq2) + 1):
        for j in range(len(seq1) + 1):
            if(i == 0):  # first row, each letter in seq1 with "-"
                if(j != 0):
                    letter_col = seq1[j - 1]
                    cur_score = result_table[i][j - 1] + int(score_dictionary.get((letter_col, "-")))
                    if type == "global":
                        result_table[i][j] = cur_score
                    elif type == "local":
                        result_table[i][j] = max(0, cur_score)
                        if result_table[i][j] > max_cell[0]:
                            max_cell = (max(max_cell[0], result_table[i][j]), i, j)
            elif(i != 0):
                letter_col = seq1[j - 1]
                letter_row = seq2[i - 1]
                if(j != 0):
                    up = result_table[i - 1][j] + int(score_dictionary.get(("-", letter_row)))
                    left = result_table[i][j - 1] + int(score_dictionary.get((letter_col, "-")))
                    up_and_left = result_table[i - 1][j - 1] + int(score_dictionary.get((letter_row, letter_col)))
                    cur_score = max(up, left, up_and_left)
                    if type == "global":
                        result_table[i][j] = cur_score
                    elif type == "local":
                        result_table[i][j] = max(0, cur_score)
                        if result_table[i][j] > max_cell[0]:
                            max_cell = (max(max_cell[0], result_table[i][j]), i, j)
                    came_from = max(up, left, up_and_left)
                    if type == "local":
                        came_from = max(came_from, 0)
                    if(came_from == up_and_left):
                        map_matrix[i][j] = (letter_col, letter_row, "up_and_left")
                    elif(came_from == up):
                        map_matrix[i][j] = ("-", letter_row, "up")
                    else:
                        map_matrix[i][j] = (letter_col, "-", "left")
                elif(j == 0):  # first column, each letter in seq2 with "-"
                    cur_score = result_table[i - 1][j] + int(score_dictionary.get(("-", letter_col)))
                    if type == "local":
                        result_table[i][j] = max(cur_score, 0)
                        if result_table[i][j] > max_cell[0]:
                            max_cell = (max(max_cell[0], result_table[i][j]), i, j)

                    else:
                        result_table[i][j] = cur_score

    if type == "global":
        return result_table[-1][-1], map_matrix
        # return result_table, map_matrix

    else:
        return max_cell, map_matrix, result_table


def global_type(fasta1, fasta2, _score_matrix):

    seq1 = read_fasta_file(fasta1)[0]       # Is this[0] because we might take more than one seq?
    seq2 = read_fasta_file(fasta2)[0]
    score_matrix = read_tsv_file(_score_matrix)
    score_dictionary = get_score_dictionary(score_matrix)
    final_score, map_matrix = calculate_score_and_path(seq1, seq2, score_dictionary, "global")
    # ----------------------------------------------------------------------------------------------
    final_score_matanya, map_matrix_matanya = calculate_score_and_path_matanya(seq1, seq2, score_dictionary, "global")
    # matrix_printer(map_matrix)
    # print("\nNow Matanya")
    # matrix_printer(map_matrix_matanya)

    # ----------------------------------------------------------------------------------------------

    res_seq1, res_seq2 = path_recover(map_matrix)
    print_results(res_seq1, res_seq2, "global", final_score)

    print("\nNow Matanya")
    res_seq1_matanya, res_seq2_matanya = path_recover(map_matrix_matanya)
    print_results(res_seq1_matanya, res_seq2_matanya, "global", final_score_matanya)



def path_recover(path_table):
    res_seq1 = ""
    res_seq2 = ""
    i = len(path_table) - 1
    j = len(path_table[0]) - 1
    while i > 0 or j >0:
        res_seq1 = path_table[i][j][0] + res_seq1
        res_seq2 = path_table[i][j][1] + res_seq2
        if path_table[i][j][2] == "up_and_left":
            i -= 1
            j -= 1
        elif path_table[i][j][2] == "up":
            i -= 1
        else:
            j -= 1
    return res_seq1, res_seq2


def local_type(fasta1, fasta2, _score_matrix):

    seq1 = read_fasta_file(fasta1)[0]
    seq2 = read_fasta_file(fasta2)[0]
    score_matrix = read_tsv_file(_score_matrix)
    score_dictionary = get_score_dictionary(score_matrix)
    max_cell, map_matrix, result_table =  calculate_score_and_path(seq1, seq2, score_dictionary, "local")

    res_seq1 = ""
    res_seq2 = ""
    i = max_cell[1]
    j = max_cell[2]
    while ((i > 0 or j >0) and result_table[i][j] != 0):
        res_seq1 = map_matrix[i][j][0] + res_seq1
        res_seq2 = map_matrix[i][j][1] + res_seq2
        if(map_matrix[i][j][2] == "up_and_left"):
            i -= 1
            j -= 1
        elif(map_matrix[i][j][2] == "up"):
            i -= 1
        else:
            j -= 1

    print_results(res_seq1, res_seq2, "local", max_cell[0])


def matrix_printer(mat):
    for row in mat:
        print("[ ", end='')
        print(row, end='')
        print(" ]")
    print()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    command_args = parser.parse_args()
    if command_args.align_type == 'global':
        global_type(command_args.seq_a, command_args.seq_b, command_args.score)
    elif command_args.align_type == 'local':
        local_type(command_args.seq_a, command_args.seq_b, command_args.score)

    elif command_args.align_type == 'overlap':
        raise NotImplementedError
    elif command_args.align_type == 'global_lin':
        raise NotImplementedError
    # print the best alignment and score



if __name__ == '__main__':
    main()
