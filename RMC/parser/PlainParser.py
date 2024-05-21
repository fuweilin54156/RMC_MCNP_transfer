# -*- coding:utf-8 -*-
# author: Kaiwen Li
# date: 2019-11-16
import re

import re
import numpy as np
from RMC.modifier.formatter.PlainFormatter import PlainFormatter
from RMC.model.input.refuelling import *
from RMC.model.input.Surface import *
from RMC.model.input.Geometry import *
from RMC.model.input.MacroBody import *
from RMC.model.input.Include import IncludeMaterial
from RMC.model.input.Criticality import Criticality
from RMC.model.input.CriticalitySearch import CriticalitySearch
from RMC.model.input.Burnup import Burnup
from RMC.model.input.Print import Print
from RMC.model.input.Material import *
from RMC.model.input.Mesh import *
from RMC.model.input.Tally import *
from RMC.model.input.Ptrac import *
from RMC.model.input.Binaryout import *
from RMC.model.input.ExternalSource import *
from RMC.model.input.FixedSource import FixedSource
from RMC.model.input.Physics import Physics
from RMC.model.input.base import Model as InputModel


class PlainParser:
    def __init__(self, inp):
        self.inp = inp
        # block list
        self.content = ""
        self.parsed_model = InputModel()
        self._is_parsed = False

    def _read_in(self):
        converter = PlainFormatter(self.inp)
        self.content = converter.format_to_cards()
        converter.clear()

    def _prepare(self):
        self.parsed_model.model['geometry'] = Geometry()

    @property
    def parsed(self):
        if self._is_parsed:
            return self.parsed_model
        self._read_in()
        self._prepare()
        for cards in self.content:
            # split block into single lines
            # example: UNIVERSE 1\nCell 1 -1:1 void=1
            #       to ["UNIVERSE 1", "Cell 1 -1:1 void=1"]
            cards = cards.strip('\n')
            card_list = cards.split('\n')
            card_title = card_list[0].split(' ', 2)
            card = card_title[0].upper()
            if card =='EXTERNALSOURCE':
                print("#####")

            if card == 'INCLUDE':
                if (len(card_title) == 2) \
                        and (card_title[1].upper().startswith('MATERIAL')):
                    self.parsed_model.model['includematerial'] = \
                        IncludeMaterial(card_title[1])
                else:
                    self.parsed_model.model['unparsed'].append(cards)
                continue
            # currently only UNIVERSE block supports >= 2
            if len(card_title) >= 2 and card != 'PTRAC':
                index = int(float(card_title[1]))
                options = ''
                # UNIVERSE with 'lattice', 'pitch', etc. options
                if len(card_title) >= 3:
                    options = card_title[2]

                if card == 'UNIVERSE':
                    univ = self._parse_universe(card_list[1:], int(index), options)
                    self.parsed_model.model['geometry'].add_universe(univ)
                else:
                    self.parsed_model.model['unparsed'].append(cards)
                    # raise ValueError('%s card with index can not be recognized!' % card)
            else:
                if card == 'REFUELLING':
                    refuel = self.__parse_refuelling(card_list[1:])
                    self.parsed_model.model['refuelling'] = refuel
                elif card == 'MATERIAL':
                    materialmodel = self.__parse_material(card_list[1:])
                    self.parsed_model.model['material'] = materialmodel
                elif card == 'CRITICALITY':
                    criticality = self.__parse_criticality(card_list[1:])
                    self.parsed_model.model['criticality'] = criticality
                elif card == 'BINARYOUT':
                    binaryout = self.__parse_binaryout(card_list[1:])
                    self.parsed_model.model['binaryout'] = binaryout
                elif card == 'EXTERNALSOURCE':
                    externalsourcemodel = self.__parse_externalsource(card_list[1:])
                    self.parsed_model.model['externalsource'] = externalsourcemodel
                elif card == 'PTRAC':
                    ptracmodel = self.__parse_ptrac(card_list)
                    self.parsed_model.model['ptrac'] = ptracmodel
                elif card == 'BURNUP':
                    burnup = self.__parse_burnup(card_list[1:])
                    self.parsed_model.model['burnup'] = burnup
                elif card == 'CRITICALITYSEARCH':
                    criticalitysearch = self.__parse_criticalitysearch(card_list[1:])
                    self.parsed_model.model['criticalitysearch'] = criticalitysearch
                elif card == 'TALLY':
                    tallymodel = self.__parse_tally(card_list[1:])
                    self.parsed_model.model['tally'] = tallymodel
                elif card == 'PRINT':
                    printmodel = self.__parse_print(card_list[1:])
                    self.parsed_model.model['print'] = printmodel
                elif card == 'MESH':
                    meshmodel = self.__parse_mesh(card_list[1:])
                    self.parsed_model.model['mesh'] = meshmodel
                elif card == 'SURFACE':
                    surfmodel = self.__parse_surface(card_list[1:])
                    self.parsed_model.model['surface'] = surfmodel
                elif card == 'MACROBODY':
                    bodymodel = self.__parse_body(card_list[1:])
                    self.parsed_model.model['macrobody'] = bodymodel
                elif card == 'FIXEDSOURCE':
                    fixedsource = self.__parse_fixed_source(card_list[1:])
                    self.parsed_model.model['fixedsource'] = fixedsource
                elif card == 'PHYSICS':
                    physics = self.__parse_physics(card_list[1:])
                    self.parsed_model.model['physics'] = physics
                else:
                    self.parsed_model.model['unparsed'].append(cards)
                    # raise ValueError('%s card can not be recognized!' % card_list[0])

        print("cards done")
        try:
            self.parsed_model.postprocess()
        except:
            print("Error self.parsed_model.postprocess")
            
        self._is_parsed = True
        return self.parsed_model

    @staticmethod
    def _parse_universe(content, index, options):
        univ = Universe(number=index)
        # if there are some options after the universe index.
        if options:
            options_val = PlainParser._parse_options(options, Universe.card_option_types)
            for opt in Universe.card_option_types.keys():
                if opt not in options_val:
                    options_val[opt] = None
            univ.add_options(options_val)

        # if there are cells in the universe.
        cells = []
        if len(content) > 0:
            for line in content:
                cells.append(PlainParser._parse_cell(line))
        univ.add_cells(cells)
        return univ

    @staticmethod
    def __parse_ptrac(content):
        content = content[0].split(' ', 1)[1]  # 把开头的PTRAC分离
        ptrac_options = PlainParser._parse_options(content, Ptrac.card_option_types)
        ptracmodel = Ptrac()
        ptracmodel.add_options(ptrac_options)
        return ptracmodel

    @staticmethod
    def __parse_binaryout(content):
        option = content[0].split(' ', 1)[1]
        opt_list = PlainParser._parse_options(option, Binaryout.card_option_types)
        binaryout_model = Binaryout()
        binaryout_model.add_options(opt_list)
        return binaryout_model

    @staticmethod
    def __parse_externalsource(content):
        source_list = []
        surfsrcread_list = []
        distributions = []
        h5source = None
        unparsed = ''
        for option in content:
            if option.split()[0].upper() == 'SOURCE':
                opt_lst = PlainParser._parse_options(option, Source.card_option_types)
                source = Source()
                source.add_options(opt_lst)
                source_list.append(source)
            elif option.split()[0].upper() == 'SURFSRCREAD':
                opt_lst = PlainParser._parse_options(option.split(' ', 1)[1], Surfsrcread.card_option_types)
                surfsrcread = Surfsrcread()
                surfsrcread.add_options(opt_lst)
                surfsrcread_list.append(surfsrcread)
            elif option.split()[0].upper() == 'DISTRIBUTION':
                distribution = PlainParser.__parse_distribution(option)
                distributions.append(distribution)
            elif option.split()[0].upper() == 'H5SOURCEINFO':
                h5source = PlainParser._parse_options(option.split(' ', 1)[1],
                                                      ExternalSource.card_option_types['H5SOURCEINFO'])
            else:
                unparsed += option + '\n'
        return ExternalSource(source=source_list, surfsrcread=surfsrcread_list, distributions=distributions,
                              h5source=h5source,
                              unparsed=unparsed)

    @staticmethod
    def __parse_distribution(content):
        id = 0
        depend = None
        type = None
        value = None
        probability = None
        bias = None

        options = content.split(' ', 2)
        if not options[0].upper() == 'DISTRIBUTION':
            raise ValueError('CELL options not found in %s' % content)

        id = int(float(options[1]))
        options = options[2].strip()

        option_parsed = PlainParser._parse_options(options, Distribution.card_option_types)
        if 'DEPEND' in option_parsed:
            depend = option_parsed['DEPEND']
        if 'TYPE' not in option_parsed:
            raise TypeError('TYPE need to be provided in Distribution card!\n')
        type = option_parsed['TYPE']
        if 'VALUE' in option_parsed:
            value = option_parsed['VALUE']
        if 'PROBABILITY' in option_parsed:
            probability = option_parsed['PROBABILITY']
        if 'BIAS' in option_parsed:
            bias = option_parsed['BIAS']
        return Distribution(id=id, depend=depend, type=type, value=value, probability=probability, bias=bias)

    @staticmethod
    def _parse_cell(content):
        options = content.split(' ', 2)
        if not options[0].upper() == 'CELL':
            raise ValueError('CELL options not found in %s' % content)
        # Get the index of the cell.
        index = int(float(options[1]))
        cell = Cell(number=index)

        # Delete the parsed part of the content.
        options = options[2].strip()

        # Get the bounds of the cell.
        cell_bounds_pattern = re.compile(r'^[^A-Za-z]+')
        m = cell_bounds_pattern.search(options)
        if not m:
            raise ValueError('Position definition not found in cell definition %s' % content)
        bounds = m.group(0)
        cell.add_bounds(bounds.strip())

        options = re.sub(r'IMP : (?P<char>[^ ])', 'IMP:\g<char>', options.upper())

        # Parse the remaining options of the cell.
        options_val = PlainParser._parse_options(options[len(bounds):], Cell.card_option_types)
        for opt in Cell.card_option_types.keys():
            if opt not in options_val:
                options_val[opt] = None
        cell.add_options(options_val)

        return cell

    @staticmethod
    def __parse_refuelling(content):
        file_option = None
        step_list = []
        index_list = []
        for option in content:
            if option.split()[0].upper() == 'FILE':
                file_option = option
            elif option.split()[0].upper() == 'REFUEL':
                options = PlainParser._parse_options(option.split(' ', 1)[1], RefuelBlock.card_option_types)
                step_list = options['STEP']
                index_list = options['INDEX']
        if file_option is None or step_list is []:
            raise ValueError('FILE option not found in refuelling block.')
        opt_list = file_option.split(' ', 1)
        return RefuelBlock(file_name=opt_list[1].strip(), steps=step_list, indexes=index_list)

    @staticmethod
    def __parse_burnup(content):
        burncell = None
        activation_cell = None
        time_step = []
        burnup_step = []
        power = []
        powerden = []
        source_rate = []
        substep = None
        inherent = None
        acelib = None
        strategy = None
        solver = None
        parallel = None
        deplete_mode = None
        xeequilibrium = None
        outputcell = None
        varymat = []
        varysurf = []
        merge = None
        succession = None
        naecf = None
        impnuc = None
        fixednuc = None
        reactions = []
        decay_source = None
        unparsed = ''
        for option in content:

            if option.split()[0].upper() == 'BURNCELL':
                burncell = PlainParser._parse_options(option, Burnup.card_option_types)['BURNCELL']
            elif option.split()[0].upper() == 'ACTIVATIONCELL':
                activation_cell = PlainParser._parse_options(option, Burnup.card_option_types)['ACTIVATIONCELL']
            elif option.split()[0].upper() == 'TIMESTEP':
                # option 不用 split是因为TIMESTEP卡中没有 m = n 这种形式
                time_step = PlainParser._parse_options(option, Burnup.card_option_types)['TIMESTEP']
            elif option.split()[0].upper() == 'BURNUPSTEP':
                burnup_step = PlainParser._parse_options(option, Burnup.card_option_types)['BURNUPSTEP']
            elif option.split()[0].upper() == 'POWER':
                power = PlainParser._parse_options(option, Burnup.card_option_types)['POWER']
            elif option.split()[0].upper() == 'POWERDEN':
                powerden = PlainParser._parse_options(option, Burnup.card_option_types)['POWERDEN']
            elif option.split()[0].upper() == 'SOURCEINTENSITY':
                source_rate = PlainParser._parse_options(option, Burnup.card_option_types)['SOURCEINTENSITY']
            elif option.split()[0].upper() == 'SUBSTEP':
                substep = PlainParser._parse_options(option, Burnup.card_option_types)['SUBSTEP']
            elif option.split()[0].upper() == 'INHERENT':
                inherent = PlainParser._parse_options(option, Burnup.card_option_types)['INHERENT']
            elif option.split()[0].upper() == 'ACELIB':
                acelib = PlainParser._parse_options(option, Burnup.card_option_types)['ACELIB']
            elif option.split()[0].upper() == 'STRATEGY':
                strategy = PlainParser._parse_options(option.split(' ', 1)[1], Burnup.card_option_types['STRATEGY'])
            elif option.split()[0].upper() == 'SOLVER':
                solver = PlainParser._parse_options(option, Burnup.card_option_types)['SOLVER']
            elif option.split()[0].upper() == 'PARALLEL':
                parallel = PlainParser._parse_options(option, Burnup.card_option_types)['PARALLEL']
            elif option.split()[0].upper() == 'DEPLETEMODE':
                deplete_mode = PlainParser._parse_options(option, Burnup.card_option_types)['DEPLETEMODE']
            elif option.split()[0].upper() == 'XEEQUILIBRIUM':
                xeequilibrium = PlainParser._parse_options(option.split(' ', 1)[1],
                                                           Burnup.card_option_types['XEEQUILIBRIUM'])
            elif option.split()[0].upper() == 'OUTPUTCELL':
                outputcell = PlainParser._parse_options(option, Burnup.card_option_types)['OUTPUTCELL']
            elif option.split()[0].upper() == 'MERGE':
                merge = PlainParser._parse_options(option.split(' ', 1)[1], Burnup.card_option_types['MERGE'])
            elif option.split()[0].upper() == 'VARYMAT':
                option_tmp = PlainParser._parse_dict_info(option.split(' ', 1)[1])
                for opt in option_tmp:
                    varymat.append(PlainParser._parse_options(opt, Burnup.card_option_types['VARYMAT']))
            elif option.split()[0].upper() == 'VARYSURF':
                option_tmp = PlainParser._parse_dict_info(option.split(' ', 1)[1])
                for opt in option_tmp:
                    varysurf.append(PlainParser._parse_options(opt, Burnup.card_option_types['VARYSURF']))
            elif option.split()[0].upper() == 'REACTIONS':
                option_tmp = PlainParser._parse_dict_info(option.split(' ', 1)[1], separator='NUCLIDE')
                for opt in option_tmp:
                    reactions.append(PlainParser._parse_options(opt, Burnup.card_option_types['REACTIONS']))
            elif option.split()[0].upper() == 'DECAYSOURCE':
                decay_source = PlainParser._parse_options(option.split(' ', 1)[1],
                                                          Burnup.card_option_types['DECAYSOURCE'])
            elif option.split()[0].upper() == 'SUCCESSION':
                succession = PlainParser._parse_options(option.split(' ', 1)[1], Burnup.card_option_types['SUCCESSION'])
            elif option.split()[0].upper() == 'NAECF':
                naecf = PlainParser._parse_options(option, Burnup.card_option_types)['NAECF']
            elif option.split()[0].upper() == 'IMPNUC':
                impnuc = PlainParser._parse_options(option, Burnup.card_option_types)['IMPNUC']
            elif option.split()[0].upper() == 'FIXEDNUC':
                fixednuc = PlainParser._parse_options(option, Burnup.card_option_types)['FIXEDNUC']
            else:
                unparsed += option + '\n'
        if burncell is []:
            raise ValueError('BURNCELL option not found in BURNUP block.')
        if time_step is [] and burnup_step is []:
            raise ValueError('TIMESTEP or BURNUPSTEP option not found in BURNUP block.')
        if time_step and burnup_step:
            raise ValueError('TIMESTEP and BURNUPSTEP are repeatly defined in BURNUP block')
        if power is [] and powerden is []:
            raise ValueError('POWER or POWERDEN option not found in BURNUP block.')
        if power and powerden:
            raise ValueError('POWER and POWERDEN are repeatly defined in BURNUP block')
        return Burnup(burn_cell=burncell, activation_cell=activation_cell, time_step=time_step, burnup_step=burnup_step,
                      power=power, power_den=powerden, source_intensity=source_rate, substep=substep, inherent=inherent,
                      acelib=acelib, strategy=strategy, solver=solver, deplete_mode=deplete_mode, parallel=parallel,
                      xeequilibrium=xeequilibrium, outputcell=outputcell, varymat=varymat, varysurf=varysurf,
                      succession=succession, merge=merge, naecf=naecf, impnuc=impnuc, fixednuc=fixednuc,
                      reactions=reactions, decay_source=decay_source,
                      unparsed=unparsed)

    @staticmethod
    def __parse_criticality(content):
        couple = None
        power_iter = None
        initsrc = None
        unparsed = ''
        for option in content:
            if option.split()[0].upper() == 'COUPLE':
                couple = PlainParser._parse_options(option.split(' ', 1)[1],
                                                    Criticality.card_option_types["COUPLE"])
            elif option.split()[0].upper() == 'POWERITER':
                power_iter = PlainParser._parse_options(option.split(' ', 1)[1],
                                                        Criticality.card_option_types["POWERITER"])
            elif option.split()[0].upper() == 'INITSRC':
                initsrc = PlainParser._parse_options(option.split(' ', 1)[1],
                                                        Criticality.card_option_types["INITSRC"])
            else:
                unparsed += option + '\n'
        return Criticality(couple=couple, power_iter=power_iter, initsrc=initsrc, unparsed=unparsed)

    @staticmethod
    def __parse_criticalitysearch(content):
        unparsed = ''
        for option in content:
            unparsed += option + '\n'
        return CriticalitySearch(unparsed=unparsed)

    @staticmethod
    def __parse_material(content):
        unparsed = ''
        mats = []
        ceace = {}
        mtlib = {}
        nubar = {}
        mgace = {}

        for option in content:
            option = re.sub(r'D (?P<char>[^ ])', 'D\g<char>', option)
            if re.search(r"MAT|SAB|DYNAMICMAT", option.split()[0].upper()):
                mats.append(PlainParser.__parse_single_mat(option))
            elif option.split()[0].upper() == "CEACE":
                ceace = PlainParser._parse_options(option.split(' ', 1)[1], Materials.card_option_types['CEACE'])
            elif option.split()[0].upper() == "MTLIB":
                mtlib = PlainParser._parse_options(option.split(' ', 1)[1], Materials.card_option_types['MTLIB'])
            elif option.split()[0].upper() == "MGACE":
                mgace = PlainParser._parse_options(option.split(' ', 1)[1], Materials.card_option_types['MGACE'])
            elif option.split()[0].upper() == "NUBAR":
                nubar = PlainParser._parse_options(option.split(' ', 1)[1], Materials.card_option_types['NUBAR'])
            else:
                unparsed += option + '\n'
        # mats.sort(key=id)
        return Materials(mats=mats, ceace=ceace, mtlib=mtlib, nubar=nubar, mgace=mgace, unparsed=unparsed)

    @staticmethod
    def __parse_single_mat(mat_card):
        opts = mat_card.split()
        mat = None
        if opts[0].upper() == "MAT":
            mat = Material(mat_id=int(opts[1]), density=opts[2])

            i = 3
            while i < len(opts):
                mat.update_nuclide(Nuclide(opts[i], opts[i + 1]))
                i += 2
        elif opts[0].upper() == "SAB":
            mat = SabMaterial(mat_id=int(opts[1]), sab_nuclides=opts[2:])
        elif opts[0].upper() == "DYNAMICMAT":
            options = PlainParser._parse_options(mat_card.split(' ', 2)[-1], DynamicMaterial.card_option_types)
            mat = DynamicMaterial(mat_id=int(opts[1]), options=options)
        return mat

    @staticmethod
    def __parse_surface(content):
        """
        Parameters
        ----------
        content: list of string
            surface信息列表

        Returns
        -------
        Surface object

        """
        surfaces = []
        for option in content:
            surf_number = int(option.split(' ', 3)[1])
            surf_type = option.split(' ', 3)[2].upper()
            surf = PlainParser._parse_options(option.split(' ', 2)[2], Surface.surf_type_para)
            surf_parameter = surf[surf_type]
            if 'BC' in surf:
                surf_boundary = surf['BC']
            else:
                surf_boundary = None
            if 'PAIR' in surf:
                surf_pair = surf['PAIR']
            else:
                surf_pair = None
            if 'TIME' in surf:
                surf_time = surf['TIME']
                surf_value = surf['VALUE']
            else:
                surf_value = None
                surf_time = None
            if 'ROTATE' in surf:
                surf_rotate = surf['ROTATE']
                surf_rotate = list(np.array(surf_rotate).reshape(3, 3))
            else:
                surf_rotate = None

            if 'ROTATEANGLE' in surf:
                surf_rotate_angle = surf['ROTATEANGLE']
            else:
                surf_rotate_angle = None

            if 'MOVE' in surf:
                surf_move = surf['MOVE']
            else:
                surf_move = None
            # 添加面到surfaces中
            # change surface type 'X/X' to 'X_X' to get the materialized object using eval() function
            surface = Surface.externalization(number=surf_number, type=surf_type, parameters=surf_parameter,
                                              boundary=surf_boundary, pair=surf_pair, time=surf_time, value=surf_value,
                                              rotate=surf_rotate, rotate_angle=surf_rotate_angle, move=surf_move)
            # 旋转面
            surface = surface.transfer()
            surfaces.append(surface)
        return Surfaces(surfaces=surfaces)

    @staticmethod
    def __parse_body(content):
        macrobodies = {}
        rotate = None
        rotate_angle = None
        move = None

        for option in content:
            body_number = int(option.split(' ', 2)[1])
            body_type = option.split(' ', 3)[2].upper()
            body_dict = PlainParser._parse_options(option.split(' ', 2)[2], MacroBody.body_type_para)
            body_parameter = body_dict[body_type]
            if 'ROTATE' in body_dict.keys():
                rotate = body_dict['ROTATE']
                rotate = np.array(rotate).reshape(3, 3)
            if 'ROTATEANGLE' in body_dict.keys():
                rotate_angle = body_dict['ROTATEANGLE']
            if 'MOVE' in body_dict.keys():
                move = body_dict['MOVE']
            if body_type == 'ELL':
                macrobodies[body_number] = ELL(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'RPP':
                macrobodies[body_number] = RPP(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'RCC':
                macrobodies[body_number] = RCC(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'BOX':
                macrobodies[body_number] = BOX(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'SPH':
                macrobodies[body_number] = SPH(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'TORUS':
                macrobodies[body_number] = TORUS(number=body_number, type=body_type, params=body_parameter,
                                                 rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'HEX' or body_type == 'RHP':
                macrobodies[body_number] = HEX(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'REC':
                macrobodies[body_number] = REC(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'TRC':
                macrobodies[body_number] = TRC(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'WED':
                macrobodies[body_number] = WED(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
            elif body_type == 'SEC':
                macrobodies[body_number] = SEC(number=body_number, type=body_type, params=body_parameter,
                                               rotate=rotate, rotate_angle=rotate_angle, move=move)
        return Macrobodies(macrobodies=macrobodies)

    @staticmethod
    def __parse_mesh(content):
        meshmodel = Mesh()
        for option in content:
            mesh_options = PlainParser._parse_options(option, MeshInfo.card_option_types)
            meshinfocard = MeshInfo()
            meshinfocard.add_options(mesh_options)
            meshmodel.add_one_mesh(meshinfocard)
        return meshmodel

    @staticmethod
    def __parse_tally(content):
        unparsed = ''
        meshtally_list = []
        surftally_list = []
        celltally_list = []
        bin_list = []
        for option in content:
            if option.split()[0].upper() == 'MESHTALLY':
                opt_lst = PlainParser._parse_options(option, MeshTally.card_option_types)
                meshtally = MeshTally()
                meshtally.add_options(opt_lst)
                meshtally_list.append(meshtally)
            elif option.split()[0].upper() == 'SURFTALLY':
                opt_lst = PlainParser._parse_options(option, SurfTally.card_option_types)
                surftally = SurfTally()
                surftally.add_options(opt_lst)
                surftally_list.append(surftally)
            elif option.split()[0].upper() == 'CELLTALLY':
                # process cell, for example: 2 > 3 > 4 *12
                option = re.sub(r'\* (?P<cell>[^ ])', '*\g<cell>', option)

                opt_lst = PlainParser._parse_options(option, CellTally.card_option_types)
                celltally = CellTally()
                celltally.add_options(opt_lst)
                celltally_list.append(celltally)
            elif option.split()[0].upper() == 'BIN':
                opt_lst = PlainParser._parse_options(option, Bin.card_option_types)
                bin_model = Bin()
                bin_model.add_options(opt_lst)
                bin_list.append(bin_model)
            else:
                unparsed += option + '\n'
        return Tally(meshtally=meshtally_list, celltally=celltally_list, surftally=surftally_list, tally_bin=bin_list,
                     unparsed=unparsed)

    @staticmethod
    def __parse_fixed_source(content):
        particle = None
        cutoff = None
        load_balance = None
        unparsed = ''
        for option in content:
            if option.split()[0].upper() == 'PARTICLE':
                particle = PlainParser._parse_options(option.split(' ', 1)[1],
                                                      FixedSource.card_option_types['PARTICLE'])
            elif option.split()[0].upper() == 'CUTOFF':
                cutoff = PlainParser._parse_options(option.split(' ', 1)[1], FixedSource.card_option_types['CUTOFF'])
            elif option.split()[0].upper() == 'LOAD_BALANCE':
                load_balance = PlainParser._parse_options(option.split(' ', 1)[1],
                                                          FixedSource.card_option_types['LOAD_BALANCE'])
            else:
                unparsed += option + '\n'
        return FixedSource(particle=particle, cutoff=cutoff, load_balance=load_balance)

    @staticmethod
    def __parse_physics(content):
        particle_mode = None
        neutron = None
        photon = None
        electron = None
        unparsed = ''
        for option in content:
            if option.split()[0].upper() == 'PARTICLEMODE':
                particle_mode = PlainParser._parse_options(option, Physics.card_option_types)['PARTICLEMODE']
            elif option.split()[0].upper() == 'NEUTRON':
                neutron = PlainParser._parse_options(option.split(' ', 1)[1],
                                                     Physics.card_option_types['NEUTRON'])
            elif option.split()[0].upper() == 'PHOTON':
                photon = PlainParser._parse_options(option.split(' ', 1)[1],
                                                    Physics.card_option_types['PHOTON'])
            elif option.split()[0].upper() == 'ELECTRON':
                electron = PlainParser._parse_options(option.split(' ', 1)[1],
                                                      Physics.card_option_types['ELECTRON'])
            else:
                unparsed += option + '\n'
        return Physics(particle_mode=particle_mode, neutron=neutron, photon=photon, electron=electron)

    @staticmethod
    def __parse_print(content):
        unparsed = ''
        inpfile = 0
        for option in content:
            if option.split()[0].upper() == 'INPFILE':
                inpfile = PlainParser._parse_options(option, Print.card_option_types)['INPFILE']
            else:
                unparsed += option + '\n'
        return Print(inpfile=inpfile, unparsed=unparsed)

    #  str: step = ... mat = ... newmat ... step ... -> [step = ... mat = ... newmat = ... , step = ... ... ...]
    @staticmethod
    def _parse_dict_info(varyinfo, separator='STEP'):
        vary_info = varyinfo.upper().split(separator)
        for i in range(1, len(vary_info)):
            vary_info[i] = separator + vary_info[i]
        vary_info.pop(0)
        return vary_info

    @staticmethod
    def _parse_options(content, cards):
        options = content.replace(' = ', ' ').split()
        options_dict = {}
        while options:
            options[0] = options[0].upper()
            if options[0] in cards:
                dtype = cards[options[0]]
                if dtype[0] == 'list':
                    [opt_val, index] = PlainParser._parse_list(options, dtype[1])
                else:
                    [opt_val, index] = PlainParser._parse_val(options, dtype[0])

                options_dict[options[0]] = opt_val
                # 移除已经解析过的部分
                options = options[index:]
            else:
                raise ValueError('%s can not be recognized in %s.' % (options[0], content))
        return options_dict

    @staticmethod
    def _parse_list(options, func=str):
        index = 1
        opt_list = []
        while index < len(options):
            # 在燃耗中可以输入星号*，后面接一正整数n，表示重复n次
            if options[index] == '*':
                try:
                    value = func(options[index - 1])
                    # todo 容错处理，如用户输入 3 * 1.2，或者 * 前后不完整
                    # todo RMC中的 * 有两种用途：一种表示重复次数，还有一种用在栅元计数器中，
                    #  表示展开所有底层栅元。现在代码还没有处理第二种逻辑。
                    # range里面 -1 是因为在读到 * 之前已经append过一次了
                    for repeat in range(int(options[index + 1]) - 1):
                        opt_list.append(value)
                    index += 2
                except ValueError:
                    option = options[index] + options[index + 1]  # 处理celltally中*36这种情况
                    opt_list.append(option)
                    index += 2
                    break
            elif re.fullmatch(r'\*\d+|\(|\)|:|>', options[index]):
                opt_list.append(options[index])
                index += 1
            elif re.fullmatch(r'[BD]', options[index].upper()):
                option = options[index] + options[index + 1]
                opt_list.append(option)
                index += 2
            else:
                try:
                    value = func(options[index])
                    index += 1
                    opt_list.append(value)
                except ValueError:
                    break
        # todo: check
        return [opt_list, index]

    @staticmethod
    def _parse_val(options, func=str):
        if func == bool:
            return [func(float(options[1])), 2]
        elif func == int:  # 二院项目要求智能识别int和float
            return [func(float(options[1])), 2]
        # todo: str process
        elif func == "concat":
            return [' '.join(options[1:]), 2]
        else:
            return [func(options[1]), 2]
