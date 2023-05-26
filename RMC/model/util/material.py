#!/usr/env/bin python3
# -*- encoding=utf-8 -*-


def get_negative_mat_file_map(univ):
    mat_dict = {}
    univ_stack = [univ]
    visited = set()
    while len(univ_stack) > 0:
        universe = univ_stack.pop()
        if universe in visited:
            continue
        visited.add(universe.number)
        if len(universe.include) > 0:
            univ_stack.extend(universe.include)
        else:
            for cell in universe:
                if cell.void:
                    continue
                if cell.fill is not None:
                    univ_stack.append(cell.include)
                elif cell.material[0] < 0:
                    mat_dict[cell.number] = -cell.material[0]
    return mat_dict
