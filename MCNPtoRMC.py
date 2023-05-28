# -*- coding:utf-8 -*-
# author: Shen PF
# date: 2021-07-23


import RMC.model.input.Material as RMCMat
import RMC.model.input.Geometry as RMCGeometry
import RMC.model.input.MacroBody as RMCMacrobody
import RMC.model.input.Criticality as RMCCriticality
from RMC.model.input.base import Model as RMCModel
import MCNP.parser.PlainParser as MCNPParser
import numpy as np
import re


def transfer(inp_MCNP):
    print('Processing file: ' + inp_MCNP + ' ...')

    M_model = MCNPParser.PlainParser(inp_MCNP).parsed
    with open(inp_MCNP + '_parsed_MCNP_model', 'w+') as f:
        f.write(str(M_model))

    R_model = RMCModel()

    # transfer surface block
    R_surfaces = []
    for M_surf in M_model.model['surface'].surfaces:
        if M_surf.tr is not None:
            if M_surf.tr.rotate is not None:
                surf = RMCGeometry.Surface.externalization(number=M_surf.number, type=M_surf.type,
                                                           parameters=M_surf.parameters,
                                                           boundary=M_surf.boundary, pair=M_surf.pair,
                                                           move=M_surf.tr.move,
                                                           rotate=np.array(M_surf.tr.rotate).reshape([3, 3]))
            else:
                surf = RMCGeometry.Surface.externalization(number=M_surf.number, type=M_surf.type,
                                                           parameters=M_surf.parameters,
                                                           boundary=M_surf.boundary, pair=M_surf.pair,
                                                           move=M_surf.tr.move)
            R_surfaces.append(surf)
        else:
            surf = RMCGeometry.Surface.externalization(number=M_surf.number, type=M_surf.type,
                                                       parameters=M_surf.parameters,
                                                       boundary=M_surf.boundary, pair=M_surf.pair)
            R_surfaces.append(surf)

    R_surfaces_model = RMCGeometry.Surfaces(surfaces=R_surfaces)
    test = str(R_surfaces_model)

    # transfer macrobody block
    R_macrobodys = []
    for M_body in M_model.model['surface'].macrobodys:
        if M_body.tr is not None:
            if M_body.tr.rotate is not None:
                body = RMCMacrobody.MacroBody.externalization(number=M_body.body_number, type=M_body.type,
                                                              parameters=M_body.params, move=M_body.tr.move,
                                                              rotate=np.array(M_body.tr.rotate).reshape([3, 3]))
            else:
                body = RMCMacrobody.MacroBody.externalization(number=M_body.body_number, type=M_body.type,
                                                              parameters=M_body.params,
                                                              move=M_body.tr.move)
            R_macrobodys.append(body)
        else:
            body = RMCMacrobody.MacroBody.externalization(number=M_body.body_number, type=M_body.type,
                                                          parameters=M_body.params)
            R_macrobodys.append(body)

    R_macrobodys_model = RMCMacrobody.Macrobodies(macrobodies=R_macrobodys)
    test = str(R_macrobodys_model)

    # transfer material block
    R_materials = []
    for cell in M_model.model['geometry'].cells:
        tmp_mat_id = cell.material
        if tmp_mat_id == 0:
            continue
        for index in range(len(M_model.model['materials'].mats)):
            if M_model.model['materials'].mats[index].mat_id == tmp_mat_id:
                M_model.model['materials'].mats[index].densities.append(cell.density)

    duplicate_mats = []
    for mat in M_model.model['materials'].mats:
        if mat.densities and len(set(mat.densities)) == 1:
            R_mat = RMCMat.Material(mat_id=mat.mat_id, density=mat.densities[0], nuclides=mat.nuclides)
            R_materials.append(R_mat)
        elif len(set(mat.densities)) > 1:
            duplicate_mats.append(mat.mat_id)
            R_mat = RMCMat.Material(mat_id=mat.mat_id, density=mat.densities[0], nuclides=mat.nuclides)
            R_materials.append(R_mat)
    if duplicate_mats:
        print(' Warning: find duplicated mat densities, id:' + str(duplicate_mats))

    R_sabs = ''
    for mt in M_model.model['materials'].mts:
        R_sabs += 'sab ' + str(mt.id) + ' ' + mt.name + '\n'

    R_materials_model = RMCMat.Materials(mats=R_materials, unparsed=R_sabs)
    test2 = str(R_materials_model)

    # transfer geometry block
    R_cells = []
    R_universes = []
    R_universes_ids = []
    for cell in M_model.model['geometry'].cells:
        out_universe_id = 0
        if cell.universe:
            out_universe_id = cell.universe
        if cell.trcl is not None and cell.lat is None:
            R_cell = RMCGeometry.Cell(number=cell.number, bounds=cell.bounds.replace('#', '!'),
                                      material=[cell.material],
                                      fill=cell.fill, imp_n=cell.impn, imp_p=cell.impp,
                                      transformation=RMCGeometry.Transformation(move=cell.trcl.move,
                                                                                rotate=cell.trcl.rotate))
        elif cell.lat is not None:
            R_cell = RMCGeometry.Cell(number=cell.number, bounds=cell.bounds.replace('#', '!'),
                                      material=[cell.material], imp_n=cell.impn, imp_p=cell.impp)
        else:
            R_cell = RMCGeometry.Cell(number=cell.number, bounds=cell.bounds.replace('#', '!'),
                                      material=[cell.material],
                                      fill=cell.fill, imp_n=cell.impn, imp_p=cell.impp)

        established_univ = False
        for index in range(len(R_universes_ids)):
            if R_universes_ids[index] == out_universe_id:
                R_universes[index].cells.append(R_cell)
                established_univ = True

        if not established_univ:
            R_lattice = None
            if cell.lat is not None and cell.lat == 1:
                [scope, pitch, fill, move] = transfer_lat1(cell, R_surfaces_model, R_macrobodys_model)
                R_lattice = RMCGeometry.Lattice(type=cell.lat, scope=scope, pitch=pitch, fill=fill)
                R_universe1 = RMCGeometry.Universe(number=out_universe_id, lattice=R_lattice,
                                                   transformation=RMCGeometry.Transformation(move=move))
                R_universe2 = RMCGeometry.Universe(number=out_universe_id * 1000 + 1)
                R_universe2.cells.append(R_cell)
                R_universes.append(R_universe2)
                R_universes_ids.append(out_universe_id * 1000 + 1)
                R_universes.append(R_universe1)
                R_universes_ids.append(out_universe_id)
            elif cell.lat is not None and cell.lat == 2:
                print(' Warning: the lat params are uncorrected processed in MCNP cell ' + str(cell.number))
                R_lattice = RMCGeometry.Lattice(type=cell.lat)
                R_universe1 = RMCGeometry.Universe(number=out_universe_id, lattice=R_lattice)
                R_universes.append(R_universe1)
                R_universes_ids.append(out_universe_id)
            else:
                R_universe = RMCGeometry.Universe(number=out_universe_id, lattice=R_lattice)
                R_universe.cells.append(R_cell)
                R_universes.append(R_universe)
                R_universes_ids.append(out_universe_id)

    # 添加 lat 里填充的 Universe 的 move 卡
    for univ in R_universes:
        if univ.lattice is not None:
            move = np.array(univ.lattice.pitch) / 2
            for fill_univ_id in univ.lattice.fill:
                for univ2 in R_universes:
                    if univ2.number == fill_univ_id:
                        univ2.transformation = RMCGeometry.Transformation(move=move)
                        break

    # 调整 universe 顺序
    R_universes = sorted(R_universes, key=lambda x: x.number)

    R_geometry_model = RMCGeometry.Geometry(universes=R_universes)
    test3 = str(R_geometry_model)

    # combine RMC model
    R_model.model['geometry'] = R_geometry_model
    R_model.model['surface'] = R_surfaces_model
    R_model.model['macrobody'] = R_macrobodys_model
    R_model.model['material'] = R_materials_model

    # set the Criticality block
    power_iter = {"KEFF0": 1, "POPULATION": [10000, 50, 300], "BATCHNUM": 1}
    R_model.model['criticality'] = RMCCriticality.Criticality(power_iter=power_iter,
                                                              unparsed='InitSrc point = 0 0 0')
    R_model.model[
        'plot'] = 'PLOT Continue-calculation = 1\n' \
                  'PlotID 1 Type = slice Color = cell Pixels=10000 10000 Vertexes=-100 -100 0 100 100 0\n' \
                  'PlotID 2 type = slice color = cell pixels=10000 10000 vertexes=-100 0 -100 100 0 100'

    # output 2 files
    with open(inp_MCNP + '_parsed_RMC_model', 'w+') as f:
        f.write(str(R_model))

    print('file: [' + inp_MCNP + '] have been processed!')


def transfer_lat1(cell, R_surfaces, R_macrobodys):
    scope = np.zeros(3)
    pitch = np.zeros(3)
    fill = np.ones(1)
    move = np.zeros(3)
    params = cell.fill  # like '-8 : 8 -8 : 8 -8 : 8 1 1 1 1 1 1 '

    x_left = 0
    x_right = 0
    y_left = 0
    y_right = 0
    z_left = 0
    z_right = 0
    index = 0
    if ':' in params:
        x_left = params[0]
        x_right = params[2]
        index = 3
        if ':' in params[3:]:
            y_left = params[3]
            y_right = params[5]
            index = 6
            if ':' in params[6:]:
                z_left = params[6]
                z_right = params[8]
                index = 9
    scope = np.array([x_right - x_left + 1, y_right - y_left + 1, z_right - z_left + 1])
    # fill = np.array(params[index:].replace(cell.universe, cell.universe*1000+1))
    fill = np.array([i if i is not cell.universe else i*1000+1 for i in params[index:]])
    nums_in_bounds = re.findall(r'[0-9]+', cell.bounds)
    nums_in_bounds = [int(i) for i in nums_in_bounds]
    if len(nums_in_bounds) > 1:  # not Macrobody case
        if x_right - x_left + 1 > 1:
            x_left_surf_num = nums_in_bounds[0]
            x_right_surf_num = nums_in_bounds[1]
            x_left_surf = None
            x_right_surf = None
            for surf in R_surfaces.surfaces:
                if surf.number == x_left_surf_num:
                    x_left_surf = surf
                if surf.number == x_right_surf_num:
                    x_right_surf = surf
            if x_left_surf.type in ['P', 'PX', 'PY', 'PZ']:
                pitch[0] = abs(x_left_surf.parameters[-1] - x_right_surf.parameters[-1])
            else:
                print(' Warning: the lat params are uncorrected processed in MCNP cell ' + str(cell.number))

        if y_right - y_left + 1 > 1:
            y_left_surf_num = nums_in_bounds[2]
            y_right_surf_num = nums_in_bounds[3]
            for surf in R_surfaces.surfaces:
                if surf.number == y_left_surf_num:
                    y_left_surf = surf
                if surf.number == y_right_surf_num:
                    y_right_surf = surf
            if y_left_surf.type in ['P', 'PX', 'PY', 'PZ']:
                pitch[1] = abs(y_left_surf.parameters[-1] - y_right_surf.parameters[-1])
            else:
                print(' Warning: the lat params are uncorrected processed in MCNP cell ' + str(cell.number))

        if z_right - z_left + 1 > 1:
            z_left_surf_num = nums_in_bounds[4]
            z_right_surf_num = nums_in_bounds[5]
            for surf in R_surfaces.surfaces:
                if surf.number == z_left_surf_num:
                    z_left_surf = surf
                if surf.number == z_right_surf_num:
                    z_right_surf = surf
            if z_left_surf.type in ['P', 'PX', 'PY', 'PZ']:
                pitch[2] = abs(z_left_surf.parameters[-1] - z_right_surf.parameters[-1])
            else:
                print(' Warning: the lat params are uncorrected processed in MCNP cell ' + str(cell.number))
    else:
        # Macrobody case, only support "RPP"
        body_num = nums_in_bounds[0]
        for body in R_macrobodys.macrobodies:
            if body_num == body.body_number:
                bound_body = body
                break
        if bound_body.type in ['RPP']:
            pitch[0] = abs(bound_body.params[0]-bound_body.params[1])
            pitch[1] = abs(bound_body.params[2]-bound_body.params[3])
            pitch[2] = abs(bound_body.params[4]-bound_body.params[5])
        else:
            print(' Warning: the lat params are uncorrected processed in MCNP cell ' + str(cell.number))

    # process move
    move[0] = pitch[0] * (x_left - 0.5)
    move[1] = pitch[1] * (y_left - 0.5)
    move[2] = pitch[2] * (z_left - 0.5)

    return scope, pitch, fill, move
