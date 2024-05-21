# -*- coding:utf-8 -*-
# author: Shen PF
# date: 2021-07-23


import RMC.model.input.Material as RMCMat
import RMC.model.input.Geometry as RMCGeometry
import RMC.model.input.MacroBody as RMCMacrobody
import RMC.model.input.Criticality as RMCCriticality
from RMC.model.input.ExternalSource import *
from RMC.model.input.base import Model as RMCModel
from RMC.parser.PlainParser import PlainParser as RMCPlainParser
import MCNP.parser.PlainParser as MCNPParser
import numpy as np
import re


def transfer(inp_MCNP):
    print('Processing file: ' + inp_MCNP + ' ...')

    M_model = MCNPParser.PlainParser(inp_MCNP).parsed
    with open(inp_MCNP + '.mcnp', 'w+') as f:
        f.write(str(M_model))

    R_model = RMCModel()

    # transfer surface block
    try:
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
        if M_model.model['surface'].unparsed is not None and M_model.model['surface'].unparsed != '':
            print(" Warning: [Surface block] unparsed items: \n" + M_model.model['surface'].unparsed)
        test = str(R_surfaces_model)
    except:
        print(" Error: failed to transfer the MCNP surf block to RMC surf.")

    # transfer macrobody block
    try:
        R_macrobodys = []
        for M_body in M_model.model['surface'].macrobodys:
            if M_body.unparsed is not None and M_body.unparsed != '':
                print(" Warning: unparsed items in MCNP macrobody " + M_body.body_number + " : " + M_body.unparsed)
            if M_body.tr is not None:
                if M_body.tr.rotate is not None:
                        body = RMCMacrobody.MacroBody.externalization(number=M_body.body_number, type=M_body.type,
                                                                    parameters=M_body.params, move=M_body.tr.move,
                                                                    rotate=np.array(M_body.tr.rotate).reshape([3, 3]))
                else:
                    body = RMCMacrobody.MacroBody.externalization(number=M_body.body_number, type=M_body.type,
                                                                  parameters=M_body.params, move=M_body.tr.move)
                R_macrobodys.append(body)  
            else:
                body = RMCMacrobody.MacroBody.externalization(number=M_body.body_number, type=M_body.type,
                                                              parameters=M_body.params)
                R_macrobodys.append(body)
                
        R_macrobodys_model = RMCMacrobody.Macrobodies(macrobodies=R_macrobodys)
        test = str(R_macrobodys_model)
    except:
        print(" Error: failed to transfer the MCNP Macrobody block to RMC Macrobody.")

    # transfer material block
    try:
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
            elif not mat.densities:
                print(' Warning: no density defined in mat: ' + str(mat.mat_id)+
                      '\nPlease check if the material has not been used. \nThe density of unused materials will be set to 0.\n')
                R_mat = RMCMat.Material(mat_id=mat.mat_id, density=888, nuclides=mat.nuclides)
                R_materials.append(R_mat)
        if duplicate_mats:
            print(' Warning: find duplicated mat densities, id: ' + str(duplicate_mats)+
                  '\nPlease use various material cards to describe the same material with different densities.\n')
        if M_model.model['materials']._unparsed is not None and M_model.model['materials']._unparsed != '':
            print(" Warning: unparsed items in MCNP materials : " + M_model.model['materials']._unparsed)

        R_sabs = ''
        for mt in M_model.model['materials'].mts:
            R_sabs += 'sab ' + str(mt.id) + ' ' + mt.name + '\n'

        R_materials_model = RMCMat.Materials(mats=R_materials, unparsed=R_sabs)
        test2 = str(R_materials_model)
    except:
        print(" Error: failed to transfer the MCNP Material block to RMC Material.")

    # transfer geometry block
    try:
        R_cells = []
        R_universes = []
        R_universes_ids = []

        Max_universe_id = 1
        for cell in M_model.model['geometry'].cells:
            if cell.universe:
                # out_universe_id 栅元U号
                Max_universe_id = max(cell.universe,Max_universe_id)
        # print(Max_universe_id)

        for cell in M_model.model['geometry'].cells:
            if cell.unparsed is not None and cell.unparsed != '':
                print(" Warning: unparsed items in MCNP cell " + str(cell.number) + " : " + cell.unparsed)
            out_universe_id = 0
            if cell.impn==0:
                test=1
            if cell.universe:
                # out_universe_id 栅元U号
                out_universe_id = cell.universe 
            if cell.trcl is not None and cell.lat is None:
                R_cell = RMCGeometry.Cell(number=cell.number, bounds=cell.bounds.replace('#', '!'),
                                          material=[cell.material],
                                          fill=cell.fill, imp_n=cell.impn, imp_p=cell.impp,
                                          transformation=RMCGeometry.Transformation(move=cell.trcl.move,rotate=cell.trcl.rotate),
                                          volume=cell.vol,temperature=cell.tmp,void=cell.void,imp_e=cell.impe)
            elif cell.lat is not None:
                R_cell = RMCGeometry.Cell(number=cell.number, bounds=cell.bounds.replace('#', '!'),
                                          material=[cell.material],imp_n=cell.impn, imp_p=cell.impp,
                                          volume=cell.vol,temperature=cell.tmp,void=cell.void,imp_e=cell.impe)
            else:

                R_cell = RMCGeometry.Cell(number=cell.number, bounds=cell.bounds.replace('#', '!'),
                                          material=[cell.material],
                                          fill=cell.fill, imp_n=cell.impn, imp_p=cell.impp,
                                          volume=cell.vol,temperature=cell.tmp,void=cell.void,imp_e=cell.impe)

            established_univ = False

            #检查是否已创建UNIVERSE
            for index in range(len(R_universes_ids)):
                if R_universes_ids[index] == out_universe_id:
                    R_universes[index].cells.append(R_cell)
                    established_univ = True
                
            #如果未创建UNIVERSE，新建一个UNIVERSE
            if not established_univ:
                R_lattice = None
                if cell.lat is not None and cell.lat == 1:

                    new_universe_id=out_universe_id + Max_universe_id

                    [scope, pitch, fill, move] = transfer_lat1(cell, R_surfaces_model, R_macrobodys_model,new_universe_id)
                    R_lattice = RMCGeometry.Lattice(type=cell.lat, scope=scope, pitch=pitch, fill=fill)
                    R_universe1 = RMCGeometry.Universe(number=out_universe_id, lattice=R_lattice,
                                                       transformation=RMCGeometry.Transformation(move=move))
                
                    R_universe2 = RMCGeometry.Universe(number=new_universe_id)
                    R_universe2.cells.append(R_cell)
                    R_universes.append(R_universe2)
                    R_universes_ids.append(new_universe_id)
                    R_universes.append(R_universe1)
                    R_universes_ids.append(out_universe_id)
                elif cell.lat is not None and cell.lat == 2:

                    new_universe_id=out_universe_id + Max_universe_id

                    [scope, sita, pitch, fill, move] = transfer_lat2(cell, R_surfaces_model, R_macrobodys_model,new_universe_id)
                    R_lattice = RMCGeometry.Lattice(type=cell.lat, scope=scope, pitch=pitch, fill=fill, theta=sita)
                    R_universe1 = RMCGeometry.Universe(number=out_universe_id, lattice=R_lattice,
                                                       transformation=RMCGeometry.Transformation(move=move))
                    R_universe2 = RMCGeometry.Universe(number=new_universe_id)
                    R_universe2.cells.append(R_cell)
                    R_universes.append(R_universe2)
                    R_universes_ids.append(new_universe_id)
                    R_universes.append(R_universe1)
                    R_universes_ids.append(out_universe_id)
                # 通常情况
                else:
                    R_universe = RMCGeometry.Universe(number=out_universe_id, lattice=R_lattice)
                    R_universe.cells.append(R_cell)
                    R_universes.append(R_universe)
                    R_universes_ids.append(out_universe_id)


        # 添加 lat 里填充的 Universe 的 move 卡
        for univ in R_universes:
            if univ.lattice is not None:
                if univ.lattice.type == 1:
                    move = np.array(univ.lattice.pitch) / 2
                elif univ.lattice.type == 2:
                    move = [0, 0, 0]
                for fill_univ_id in univ.lattice.fill:
                    for univ2 in R_universes:
                        if univ2.number == fill_univ_id:
                            if univ2.transformation is not None:
                                if not (univ2.transformation.move == move).all():
                                    print(" Warning: the variables of 'move' in Universe " + str(univ2.number) + " are uncorrected processed!")
                            univ2.transformation = RMCGeometry.Transformation(move=move)
                            break

        # 调整 universe 顺序
        R_universes = sorted(R_universes, key=lambda x: x.number)
        R_geometry_model = RMCGeometry.Geometry(universes=R_universes)
        test3 = str(R_geometry_model)
        # print(test3)
    except:
        print(" Error: failed to transfer the MCNP Geometry(Cell) block to RMC Geometry.")

    # transfer the SDEF block
    try:
        if M_model.model['externalsource'] is not None:
            Combined_distribution = Combine_distrtibution(M_model.model['externalsource'].distributions)
            h5source = None
            surfsrcread_list = []
            unparsed = ''
            [RMC_Distribution, RMC_source] = Transfer_source(Combined_distribution, M_model.model['externalsource'].source, R_universes,M_model.model['mode'].mode)
            R_externalsource = ExternalSource(source=RMC_source, surfsrcread=surfsrcread_list, distributions=RMC_Distribution,
                                              h5source=h5source, unparsed=unparsed)
            R_model.model['externalsource'] = R_externalsource
    except:
        print(" Error: failed to transfer the MCNP SDEF block to RMC Source Definition.")

    # transfer the kcode block
    try:
        if M_model.model['critical'] is not None:
            power_iter = None
            initsrc = None
            if len(M_model.model['critical'].kcode) > 0:
                kcode_option = M_model.model['critical'].kcode['KCODE']
                power_iter = {"KEFF0": kcode_option[1], "POPULATION": [kcode_option[0], kcode_option[2], kcode_option[3]], "BATCHNUM": 1}
            if len(M_model.model['critical'].ksrc) > 0:
                ksrc_option = M_model.model['critical'].ksrc['KSRC']
                initsrc = {"POINT": ksrc_option}
            if M_model.model['externalsource'] is not None:
                initsrc = {"EXTERNALSOURCE": RMC_source[0]._id}
            R_model.model['criticality'] = RMCCriticality.Criticality(power_iter=power_iter, initsrc=initsrc)
    except:
        print(" Error: failed to transfer the MCNP Kcode block to RMC Criticality.")

    # combine RMC model
    try:
        R_model.model['geometry'] = R_geometry_model
        if len(M_model.model['surface'].surfaces) > 0:
            R_model.model['surface'] = R_surfaces_model
        if len(M_model.model['surface'].macrobodys) > 0:
            R_model.model['macrobody'] = R_macrobodys_model
        R_model.model['material'] = R_materials_model

        # output the Unparsed blocks
        if M_model.model['unparsed'] is not None and len(M_model.model['unparsed']) > 0:
            print(" Warning: unparsed blocks in MCNP file : " + '  '.join(M_model.model['unparsed']))

        # set the Criticality block
        # power_iter = {"KEFF0": 1, "POPULATION": [10000, 50, 300], "BATCHNUM": 1}
        # R_model.model['criticality'] = RMCCriticality.Criticality(power_iter=power_iter,
        #                                                           unparsed='InitSrc point = 0 0 0')
        # set the plot block
        R_model.model['plot'] = 'PLOT\n' \
                      'PlotID 1 Type = slice Color = cell Pixels=10000 10000 Vertexes=-100 -100 0 100 100 0\n' \
                      'PlotID 2 type = slice color = cell pixels=10000 10000 vertexes=-100 0 -100 100 0 100'

        # output 2 files, RMC python file and RMC binary file
        

        with open(inp_MCNP + '.rmc.python', 'w+') as f: 
            s=str(R_model)
            f.write(s)
        with open(inp_MCNP + '.rmc.binary', 'w+') as f:
            R_py_file = inp_MCNP + '.rmc.python'
            try:
                R_parserd_model = RMCPlainParser(R_py_file).parsed
                f.write(str(R_parserd_model))

            except:  # ValueError as e:
                print(" Warning: could not generat RMC binary input file.")
                # print(" Catch the exception: ", e)
    except:
        print(" Error: failed to output the RMC Model.")

    print('File: [' + inp_MCNP + '] have been processed!\n')


def transfer_lat1(cell, R_surfaces, R_macrobodys,new_universe_id):
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
    fill = np.array([i if i is not cell.universe else new_universe_id for i in params[index:]])
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

    # if the scope is 1, set zero to the pitch and move values
    for i in range(3):
        if scope[i] == 1:
            move[i] = 0
            pitch[i] = 0

    return scope, pitch, fill, move


def transfer_lat2(cell, R_surfaces, R_macrobodys,new_universe_id):
    scope = np.zeros(2)
    pitch = np.zeros(2)
    fill = np.ones(1)
    move = np.zeros(3)
    sita = 0
    params = cell.fill  # like '-8 : 8 -8 : 8 -8 : 8 1 1 1 1 1 1 '
    x_left = 0
    x_right = 0
    y_left = 0
    y_right = 0
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
    scope = np.array([x_right - x_left + 1, y_right - y_left + 1])
    # fill = np.array(params[index:].replace(cell.universe, cell.universe*1000+1))
    fill = np.array([i if i is not cell.universe else new_universe_id for i in params[index:]])
    nums_in_bounds = re.findall(r'[0-9]+', cell.bounds)
    nums_in_bounds = [int(i) for i in nums_in_bounds]
    if len(nums_in_bounds) > 1:  # not Macrobody case
        if x_right - x_left + 1 > 1:
            x_left_surf_num = nums_in_bounds[0]
            x_right_surf_num = nums_in_bounds[1]
            y_left_surf_num = nums_in_bounds[2]
            y_right_surf_num = nums_in_bounds[3]
            z_left_surf_num = nums_in_bounds[2]
            z_right_surf_num = nums_in_bounds[3]
            x_left_surf = None
            x_right_surf = None
            y_left_surf = None
            y_right_surf = None
            z_left_surf = None
            z_right_surf = None
            for surf in R_surfaces.surfaces:
                if surf.number == x_left_surf_num:
                    x_left_surf = surf
                if surf.number == x_right_surf_num:
                    x_right_surf = surf
                if surf.number == y_left_surf_num:
                    y_left_surf = surf
                if surf.number == y_right_surf_num:
                    y_right_surf = surf
                if surf.number == z_left_surf_num:
                    z_left_surf = surf
                if surf.number == z_right_surf_num:
                    z_right_surf = surf
            pitch[0] = abs(x_right_surf.parameters[-1] - x_left_surf.parameters[-1])
            y_length = None
            z_length = None
            if len(y_left_surf.parameters) == 4:
                y_length = abs(y_left_surf.parameters[3] - y_right_surf.parameters[3]) / math.sqrt(
                    y_left_surf.parameters[0] * y_left_surf.parameters[0] +
                    y_left_surf.parameters[1] * y_left_surf.parameters[1] +
                    y_left_surf.parameters[2] * y_left_surf.parameters[2])
            if len(z_left_surf.parameters) == 4:
                z_length = abs(z_left_surf.parameters[3] - z_right_surf.parameters[3]) / math.sqrt(
                    z_left_surf.parameters[0] * z_left_surf.parameters[0] +
                    z_left_surf.parameters[1] * z_left_surf.parameters[1] +
                    z_left_surf.parameters[2] * z_left_surf.parameters[2])
            if y_length == z_length:
                pitch[1] = y_length
            tantheta = abs(y_left_surf.parameters[0]/y_left_surf.parameters[1])
            sita = 90 - math.atan(tantheta)*180/math.pi
    else:
        # Macrobody case, only support "RPP"
        body_num = nums_in_bounds[0]
        for body in R_macrobodys.macrobodies:
            if body_num == body.body_number:
                bound_body = body
                break
        if bound_body.type in ['HEX'] or bound_body.type in ['RHP']:
            if len(bound_body.params) == 9:
                length = math.sqrt(bound_body.params[6]*bound_body.params[6]+bound_body.params[7]*bound_body.params[7]+
                                   bound_body.params[8]*bound_body.params[8])
                pitch[0] = 2 * length
                pitch[1] = 2 * length
                sita = 60
            if len(bound_body.params) == 15:
                first_vector = [bound_body.params[6], bound_body.params[7], bound_body.params[8]]
                first_length = math.sqrt(bound_body.params[6]*bound_body.params[6]+bound_body.params[7]*bound_body.params[7]+bound_body.params[8]*bound_body.params[8])
                second_vector = [bound_body.params[9], bound_body.params[10], bound_body.params[11]]
                second_length = math.sqrt(bound_body.params[9]*bound_body.params[9]+bound_body.params[10]*bound_body.params[10]+bound_body.params[11]*bound_body.params[11])
                third_vector = [bound_body.params[12], bound_body.params[13], bound_body.params[14]]
                third_length = math.sqrt(bound_body.params[12]*bound_body.params[12]+bound_body.params[13]*bound_body.params[13]+bound_body.params[14]*bound_body.params[14])
                if first_vector[1] == 0 and first_vector[2] == 0:
                    x_length = first_length
                    if second_length == third_length:
                        y_length = second_length
                        cos_theta = first_vector[0]*second_vector[0]/(x_length*y_length)
                if second_vector[1] == 0 and second_vector[2] == 0:
                    x_length = second_length
                    if first_length == third_length:
                        y_length = first_length
                        cos_theta = second_vector[0]*first_vector[0]/(x_length*y_length)
                if third_vector[1] == 0 and third_vector[2] == 0:
                    x_length = third_length
                    if first_length == second_length:
                        y_length = first_length
                        cos_theta = third_vector[0] * first_vector[0]/(x_length * y_length)
                pitch[0] = 2 * x_length
                pitch[1] = 2 * y_length
                if math.acos(cos_theta)*180/math.pi < 90:
                    sita = math.acos(cos_theta) * 180/math.pi
                else:
                    sita = 180 - math.acos(cos_theta) * 180/math.pi
        else:
            print(' Warning: the lat params are uncorrected processed in MCNP cell ' + str(cell.number))

    # process move
    move[0] = (x_left + y_left/2) * pitch[0]
    move[1] = -math.sqrt(pow(y_left * pitch[1], 2)-pow(y_left * pitch[1]/2, 2))
    move[2] = 0

    return scope, sita, pitch, fill, move


def Transfer_source(distributions, sources, R_universes,mode):

    # 转换所需结果为RMC的EXTERNALSOURCE所需5个部分
    rmc_distribution = []
    rmc_source = []
    rmc_surfsrcread = []
    rmc_transformation = []
    rmc_h5sourceinfo = []

    function_transfer = {'-2': -3, '-3': -4, '-4': -5, '-5': -6, '-21': -1, '-31': -2}


    #MCNP以 F VAR表示相关的变量，RMC以depend表示，应先遍历source，再描述distribution
    for pos, source in enumerate(sources):

        fraction_value = None # MCNP无多个源设置。如果实际使用了若干个源，需要在分布中设置合适的probability，或者最后应该手动检查概率是否归一化
        partilce_value = None

        point_value = None # RMC快速定义源设置
        sphere_value = None # RMC快速定义源设置
        cyl_x_value = None # RMC快速定义源设置
        cyl_y_value = None # RMC快速定义源设置
        cyl_z_value = None # RMC快速定义源设置

        surface_value = None
        x_value = None
        y_value = None
        z_value = None
        position_value = None
        radius_value = None
        axis_value = None
        polar_value = None
        polartheta_value = None
        extent_value = None
        height_value = None
        norm_value = None
        vector_value = None
        direcmiu_value = None
        faivector_value = None
        direcfai_value = None
        energy_value = None
        cell_value = None
        sampeff_value = None
        biasfrac_value = None
        weight_value = None
        transform_value = None

        
        # 粒子源类型
        if source._par is not None:
            partilce_value = source._par
        else :
            if mode == ['P'] :
                partilce_value = 2
            else :
                partilce_value = 1
            # partilce_value = int(1)


        # 面源参数
        if source._sur is not None:
            surface_value = source._sur
        
        # 点源参数
        if source._x is not None:
            if re.search(r'F', source._x[0]):
                x_value = [source._x[1]]
            else:
                x_value = source._x
        if source._y is not None:
            if re.search(r'F', source._y[0]):
                y_value = [source._y[1]]
            else:
                y_value = source._y
        if source._z is not None:
            if re.search(r'F', source._z[0]):
                z_value = [source._z[1]]
            else:
                z_value = source._z
        if source.pos is not None:
            position_value = source.pos

        # 体源参数
        if source._rad is not None:
            radius_value = source._rad
        if source._axs is not None:
            axis_value = source._axs
        # 面源、体源的EXT定义不同；默认为体源，EXT为柱体高度；sur非0时为面源，EXT为余弦
        if source._ext is not None:
            if source._sur:
                extent_value = source._ext
            else:
                height_value = source._ext
        
        # 粒子发射角度参数
        if source._nrm is not None:
            norm_value = source._nrm
        if source._dir is not None:
            direcmiu_value = source._dir
        if source._vec is not None:
            vector_value = source._vec

        # 粒子发射能量参数
        if source._erg is not None:
            energy_value = source._erg
        
        # 粒子所在栅元参数，注意MCNP与RMC的栅元向量写法顺序是相反的（MCNP:1<10<100 RMC:100>10>1）
        if source.cell is not None:
            cell_value = source.cell
        
        # 源粒子权重参数，只能是单个正数
        if source._wgt is not None:
            weight_value = source._wgt
        
        # 源变换参数，为trID
        if source._tr is not None:
            transform_value = source._tr

        r_Source = Source(source_id=1, fraction=fraction_value, particle=partilce_value, point=point_value, sphere=sphere_value, cyl_x=cyl_x_value, cyl_y=cyl_y_value,
                 cyl_z=cyl_z_value, surface=surface_value, x=x_value, y=y_value, z=z_value, position=position_value, radius=radius_value, axis=axis_value, polar=polar_value,
                 polartheta=polartheta_value, extent=extent_value, height=height_value, norm=norm_value, vector=vector_value, direcmiu=direcmiu_value, faivector=faivector_value,
                 direcfai=direcfai_value, energy=energy_value, cell=cell_value, sampeff=sampeff_value, biasfrac=biasfrac_value, weight=weight_value, transform=transform_value)
    
        rmc_source.append(r_Source)


    transfer_id = {}
    for source in sources:
        #RMC distribution type 0 相当于MCNP S分布
        #RMC distribution type 1 2 3相当于MCNP L分布
        #RMC distribution type 4 相当于MCNP H分布
        #RMC distribution type 5 相当于MCNP A分布
        #type 0 1 2 3 4才可用depend依赖，type5是连续分布不能被depend
        #depend的用法：var2 depend = var1，当var1取第一个值时，var2用第一种分布抽样，当var1取第二个值时，var2用第二种分布抽样

    
        #RMC distribution type 2
        if source.pos:
            Pos = source.pos
            for pos in Pos:
                # Pos为[x,y,z]的取值
                if re.match(r'D', str(pos), flags=re.IGNORECASE):
                    if 'POS' in transfer_id.keys():
                        transfer_id['POS'] += int(re.findall("\d+", pos)[0])
                    else:
                        transfer_id['POS'] = int(re.findall("\d+", pos)[0])
        #RMC distribution type 3
        if source.cell:
            cells = source.cell
            for cell in cells:
                if re.match(r'D', str(cell), flags=re.IGNORECASE):
                    if 'CELL' in transfer_id.keys():
                        transfer_id['CELL'] += int(re.findall("\d+", cell)[0])
                    else:
                        transfer_id['CELL'] = int(re.findall("\d+", cell)[0])
    


    # distribution_ID 为MCNP中的SI几，或者RMC中的distribution <ID>
    # distributions : {1: {'SI': [...], 'SP': [...]}, 2: {'SI': [...], 'SP': [...]}}
    #      RMC的distribution格式：
    #      Distribution  <id> 
    #      [Depend = <id>]  Type = <flag>  Value= <params>
    #      Probability = <params>  [Bias = <params>]

    distribution_ID = distributions.keys()
    for id in distribution_ID:
        
        type = None
        type_value = []
        probability_value = []
        bias_value = []
        depend_id = None

        cell = False
        pos = False

        if transfer_id:
            if 'CELL' in transfer_id.keys():
                cell_id = [transfer_id['CELL']]
                if id in cell_id:
                    cell = True
            if 'POS' in transfer_id.keys():
                pos_id = [transfer_id['POS']]
                if id in pos_id:
                    pos = True
        
        
        
        distribution_value = distributions[id]
        
        if "SI" in distribution_value:
            # for source in sources:
                # if source._rad is not None or source._ext is not None:
                #     if re.match(r'D\d+', source._rad[0], flags=re.I) or re.match(r'D\d+', source._ext[0], flags=re.I):
                #         dis_id = int(list(source._rad[0])[1:])
                #         if dis_id == id and "SP" not in distribution_value:
                #             type = function_transfer['-21']
                #             probability_value.append(1)
                #             type_value.append(0)

            #RMC distribution type 0 相当于MCNP S分布
            #RMC distribution type 1 2 3相当于MCNP L分布
            #RMC distribution type 4 相当于MCNP H分布
            #RMC distribution type 5 相当于MCNP A分布
            si_value = distribution_value['SI']
            if cell:
                type_value = postprocess(si_value, R_universes)
                type = 3
            si_option = None
            for value in si_value:
                if value.isalpha():
                    si_option = value
                else:
                    if not cell:
                        type_value.append(value)
            if si_option == None:
                si_option = 'H'

            if si_option == 'H':
                type = 4
            elif si_option == 'S':
                type = 0
            elif si_option == 'L':
                type = 1
            elif si_option == 'A':
                type = 5
            else:
                ValueError('Unknow Distribution Type')
        

        
        if "DS" in distribution_value:
            ds_value = distribution_value['DS']
            ds_option = None
            depend_id = None
            type = None
            pattern = r'depend_distribution_id=(\d+)'
            for item in ds_value:
                # 检查是否匹配 depend_distribution_id
                if item.startswith('depend_distribution_id='):
                    
                    match = re.search(pattern, item)
                    if match:
                        depend_id = int(match.group(1))  # 提取数字
                # 检查是否是字母
                elif item.isalpha():
                    ds_option = item  # 直接赋值
                elif item.isdigit():
                    type_value.append(int(item))

            if ds_option == None:
                ds_option = 'H'

            if ds_option == 'H':
                type = 4
            elif ds_option == 'S':
                type = 0
            elif ds_option == 'L':
                type = 1
            elif ds_option == 'A':
                type = 5
            else:
                ValueError('Unknow Distribution Type')



        if "SP" in distribution_value:
            sp_value = distribution_value['SP']
            sp_option = None
            for value in sp_value:
                if value.isalpha():
                    sp_option = value
                else:
                    if float(value) < 0:
                        type = function_transfer[value]
                    else:
                        probability_value.append(float(value))
            if sp_option == None:
                sp_option = 'D'
        if "SB" in distribution_value:
            sb_value = distribution_value['SB']
            sb_option = None
            for value in sb_value:
                if value.isalpha():
                    sb_option = value
                else:
                    bias_value.append(value)
            if sb_option == None:
                sb_option = 'D'
        # 检查 "SP" 是否不在 distribution_value 中
    if probability_value==[]:
        total_elements = len(type_value)
        
        # 避免除以零错误
        if total_elements == 0:
            raise ValueError("Source type_value list is empty, cannot divide by zero")

        # 根据 type 生成 probability_value
        if type in [0, 1, 3, 5]:  # 对于类型 0, 1, 3, 5
            probability_value = [(1 / total_elements) for _ in range(total_elements)]
        elif type == 2:  # 对于类型 2
            # 元素个数是 type_value 列表元素个数的 1/3，向下取整
            probability_value_count = total_elements // 3
            probability_value = [(1 / probability_value_count) for _ in range(probability_value_count)]
        elif type == 4:  # 对于类型 4
            # 元素个数是 type_value 列表元素个数 - 1，需要确保总数不为0
            probability_value_count = total_elements - 1 if total_elements > 0 else 0
            probability_value = [(1 / probability_value_count) for _ in range(probability_value_count)]

            

        r_Distribution=Distribution(id=id, depend=depend_id, type=type, value=type_value, probability=probability_value, bias=bias_value)
        
        
        
        rmc_distribution.append(r_Distribution)
        
    return [rmc_distribution, rmc_source]


def Combine_distrtibution(source_distribution):
    distributions = {}
    for distribution in source_distribution:
        inner_dict = {}
        if distribution.SP is not None:
            if int(distribution.SP[0]) in distributions.keys():
                inner_dict['SP'] = distribution.SP[1:len(distribution.SP)]
                distributions[int(distribution.SP[0])].update(inner_dict)
            else:
                inner_dict['SP'] = distribution.SP[1:len(distribution.SP)]
                distributions[int(distribution.SP[0])] = inner_dict
        if distribution.SI is not None:
            if int(distribution.SI[0]) in distributions.keys():
                inner_dict['SI'] = distribution.SI[1:len(distribution.SI)]
                distributions[int(distribution.SI[0])].update(inner_dict)
            else:
                inner_dict['SI'] = distribution.SI[1:len(distribution.SI)]
                distributions[int(distribution.SI[0])] = inner_dict
        if distribution.SC is not None:
            if int(distribution.SC[0]) in distributions.keys():
                inner_dict['SC'] = distribution.SC[1:len(distribution.SC)]
                distributions[int(distribution.SC[0])].update(inner_dict)
            else:
                inner_dict['SC'] = distribution.SC[1:len(distribution.SC)]
                distributions[int(distribution.SC[0])] = inner_dict
        if distribution.DS is not None:
            if int(distribution.DS[0]) in distributions.keys():
                inner_dict['DS'] = distribution.DS[1:len(distribution.DS)]
                distributions[int(distribution.DS[0])].update(inner_dict)
            else:
                inner_dict['DS'] = distribution.DS[1:len(distribution.DS)]
                distributions[int(distribution.DS[0])] = inner_dict
    return distributions


def postprocess(value, R_universes=None):
    """
    兼容MCNP风格的通用源Cell输入，即采用三点坐标表示lattice位置，如2 > (2 1 1)> > 6，并将其改为
    RMC风格的通用源输入，即使用数字序号表示lattice文职，如2 > 3 > 6 (2*2 lattice)
    """
    if not value:
        return
    value_str = ' '.join([str(int(x)) if Distribution.is_integral(x) else str(x) for x in value])

    if not re.search(r':|\(|\)', value_str):
        return
    # get all the cell expansion with (), for example:
    # expansion: 21 > 5>6>60 2 > ( 2  1 1)>( 2 3 4 ) >  6  2 > (1 1 1)>( 2 3 5)> 66
    # get [2 > ( 2  1 1)>( 2 3 4 ) >  6, 2 > (1 1 1)>( 2 3 5)> 66]
    # cell_compile = re.compile(r'\d+\s*:\s*\d+\s*\(\s*-?\d+\s*-?\d+\s*-?\d+\s*\)\s*:\s*\d+')
    cell_compile = re.compile(r'\d+\s*:(?:\s*\d+\s*\(\s*-?\d+\s*-?\d+\s*-?\d+\s*\)\s*:)*\s*\d+')
    cell_list = re.findall(cell_compile, value_str)
    # replace the expansion using the lattice index

    post_value = []
    cell_dict = {}
    univ_dict = {}
    for univ in R_universes:
        univ_dict[univ.number] = univ
        for cell in univ.cells:
            cell_dict[cell.number] = cell
    for expansion in cell_list:
        process_cell_expansion(expansion=expansion, univ_dict=univ_dict, cell_dict=cell_dict, value=post_value)
    value = post_value
    return value


def process_cell_expansion(expansion=None, univ_dict=None, cell_dict=None, value=None):
    """
    将MCNP风格的栅元展开式，修改为RMC风格的栅元展开式
    """
    if not re.search(r'\(|\)', expansion):
        return
    # split the single expansion, ['2', ' 2  1 1', ' 2 3 4 ', ' 6']
    cell_expansion = [x for x in re.split(r':|\(|\)', expansion) if x != ' ']
    scope = []
    fill = []
    while cell_expansion:
        cell_lattice = re.findall(r'-?\d+', cell_expansion[0])

        if len(cell_lattice) == 1:
            cell = cell_dict[int(cell_lattice[0])]
            # for the bottom cell without filling
            if cell.fill is None:
                if int(cell_lattice[0]) == int(re.findall(r'-?\d+', cell_expansion[-1])[0]):
                    value.append(int(cell_lattice[0]))
                del cell_expansion[0]
                continue
            else:
                # for cell with filling
                fill_universe = univ_dict[cell.fill]
                # for universe with lattice, get lattice scope
                if fill_universe.lattice is not None:
                    scope = fill_universe.lattice.scope
                    fill = fill_universe.lattice.fill

                value.append(int(cell_lattice[0]))
                value.append('>')
                del cell_expansion[0]

        elif len(cell_lattice) == 3:
            if int(cell_lattice[0]) < 0:
                x = abs(int(cell_lattice[0])) + 1
            else:
                x = 2 * int(cell_lattice[0]) + 1
            if int(cell_lattice[1]) < 0:
                y = abs(int(cell_lattice[1])) + 1
            else:
                y = 2 * int(cell_lattice[1]) + 1
            if int(cell_lattice[2]) < 0:
                z = abs(int(cell_lattice[2])) + 1
            else:
                z = 2 * int(cell_lattice[2]) + 1


            lattice_index = x + scope[0] * (y - 1) + scope[0] * scope[1] * (z - 1)

            universe_id = fill[lattice_index - 1]
            fill_universe = univ_dict[universe_id]
            if fill_universe.lattice is not None:
                scope = fill_universe.lattice.scope
                fill = fill_universe.lattice.fill

            value.append(lattice_index)
            value.append('>')
            del cell_expansion[0]



