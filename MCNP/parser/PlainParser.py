# -*- coding:utf-8 -*-
# author: ShenPF
# date: 2021-07-20

import re
import copy
from MCNP.parser.PlainFormatter import PlainFormatter
from MCNP.model.base import Model as InputModel
from MCNP.model.Geometry import *
from MCNP.model.Material import *
from MCNP.model.Source import *
from MCNP.model.Critically import *


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

    @property
    def parsed(self):
        if self._is_parsed:
            return self.parsed_model
        self._read_in()

        self.parsed_model.model['unparsed'] = []

        geometry_card = self.content[0].split('\n')
        surface_card = self.content[1].split('\n')

        try:
            [cells, geom_unparsed] = self.__parse_geometry(geometry_card)
            geometry_model = Geometry(cells=cells, unparsed=geom_unparsed)
        except:
            print(" Error: failed to parse the MCNP geometry block")

        try:
            surface_model = self.__parse_surface_macrobody(surface_card)
        except:
            print(" Error: failed to parse the MCNP Surface block")

        if len(self.content) > 1:
            other_card = self.content[2]
            other_cards = self.split_othercard(other_card)
            # process MAT card
            try:
                mat_card = other_cards[0]
                material_model = self.__parse_material(mat_card)
                self.parsed_model.model['materials'] = material_model
            except:
                print(" Error: failed to parse the MCNP Material block")
                self.parsed_model.model['unparsed'] += other_cards[0]

            # process IMP card
            try:
                imp_card = other_cards[1]
                imp_dict = self.__parse_importance(imp_card)
                if 'IMP:N' in imp_dict:
                    imp_n_list = imp_dict['IMP:N']
                    for i in range(len(geometry_model.cells)):
                        geometry_model.cells[i].impn = imp_n_list[i]
                if 'IMP:P' in imp_dict:
                    imp_p_list = imp_dict['IMP:P']
                    for i in range(len(geometry_model.cells)):
                        geometry_model.cells[i].impp = imp_p_list[i]
                if 'IMP:E' in imp_dict:
                    imp_e_list = imp_dict['IMP:E']
                    for i in range(len(geometry_model.cells)):
                        geometry_model.cells[i].impn = imp_e_list[i]
            except:
                print(" Error: failed to parse the MCNP importance definition block")
                self.parsed_model.model['unparsed'] += other_cards[1]

            # process TR card
            try:
                tr_card = other_cards[2]
                tr_model = self.__parse_tr(tr_card)
                for cell in geometry_model.cells:
                    if cell.trcl is not None:
                        if cell.trcl.num is not None:
                            for i in range(len(tr_model)):
                                if tr_model[i].num == cell.trcl.num:
                                    cell.trcl.move = tr_model[i].move
                                    cell.trcl.rotate = tr_model[i].rotate
                for surface in surface_model.surfaces:
                    if surface.tr is not None:
                        if surface.tr.num is not None:
                            for i in range(len(tr_model)):
                                if tr_model[i].num == surface.tr.num:
                                    surface.tr.move = tr_model[i].move
                                    surface.tr.rotate = tr_model[i].rotate
                for body in surface_model.macrobodys:
                    if body.tr is not None:
                        if body.tr.num is not None:
                            for i in range(len(tr_model)):
                                if tr_model[i].num == body.tr.num:
                                    body.tr.move = tr_model[i].move
                                    body.tr.rotate = tr_model[i].rotate
            except:
                print(" Error: failed to parse the MCNP TR or TRCL card")
                self.parsed_model.model['unparsed'] += other_cards[2]

            # process SDEF card
            try:
                if len(other_cards[3]) > 0:
                    sdef_card = other_cards[3]
                    sdef_model = self.__parse_sdef(sdef_card)
                    self.parsed_model.model['externalsource'] = sdef_model
            except:
                print(" Error: failed to parse the MCNP Source Definition block")
                self.parsed_model.model['unparsed'] += other_cards[3]

            # process kcode card
            try:
                if len(other_cards[4]) > 0:
                    critical_card = other_cards[4]
                    critical_model = self.__parse_critical(critical_card)
                    self.parsed_model.model['critical'] = critical_model
            except:
                print(" Error: failed to parse the MCNP kcode block")
                self.parsed_model.model['unparsed'] += other_cards[4]

            self.parsed_model.model['unparsed'] += other_cards[5]

        try:
            self.parsed_model.model['surface'] = surface_model
            self.parsed_model.model['geometry'] = geometry_model
            self.parsed_model.postprocess()
        except:
            print(" Error: failed to parse the MCNP geometry block")

        return self.parsed_model

    @staticmethod
    def __parse_importance(content):
        unparsed = ''
        imp_dict = {}
        for option in content:
            if re.match(r'IMP\:N', option, re.I):
                words = option.split()
                values = [float(words[i + 1]) for i in range(len(words) - 1)]
                imp_dict["IMP:N"] = values
            elif re.match(r'IMP\:E', option, re.I):
                words = option.split()
                values = [float(words[i + 1]) for i in range(len(words) - 1)]
                imp_dict["IMP:E"] = values
            elif re.match(r'IMP\:P', option, re.I):
                words = option.split()
                values = [float(words[i + 1]) for i in range(len(words) - 1)]
                imp_dict["IMP:P"] = values
            else:
                unparsed += option + '\n'
        return imp_dict

    @staticmethod
    def split_othercard(other_card):
        lines = other_card.split('\n')
        mat_lines = []
        imp_lines = []
        trcl_lines = []
        sdef_lines = []
        critical_lines = []
        un_parsed = []
        for line in lines:
            words = line.split(' ')
            if re.match(r'm[1-9]+', words[0], re.I):
                mat_lines.append(line)
            elif re.match(r'mt[1-9]+', words[0], re.I):
                mat_lines.append(line)
            elif re.match(r'IMP', words[0], re.I):
                imp_lines.append(line)
            elif re.match(r'tr', words[0], re.I) or re.match(r'\*tr', words[0], re.I):
                trcl_lines.append(line)
            elif re.match(r'SDEF', words[0], re.I):
                sdef_lines.append(line)
            elif re.match(r'SP', words[0], re.I):
                sdef_lines.append(line)
            elif re.match(r'SI', words[0], re.I):
                sdef_lines.append(line)
            elif re.match(r'SB', words[0], re.I):
                sdef_lines.append(line)
            elif re.match(r'DS', words[0], re.I):
                sdef_lines.append(line)
            elif re.match(r'KCODE', words[0], re.I):
                critical_lines.append(line)
            elif re.match(r'KSRC', words[0], re.I):
                critical_lines.append(line)
            else:
                un_parsed.append(line)
        return [mat_lines, imp_lines, trcl_lines, sdef_lines, critical_lines, un_parsed]

    @staticmethod
    def __parse_material(content):
        unparsed = ''
        mats = []
        mts = []
        for option in content:
            if re.match(r'm[1-9]+', option, re.I):
                mats.append(PlainParser.__parse_single_mat(option))
            elif re.match(r'mt[1-9]+', option, re.I):
                words = option.split()
                mt_id = int(words[0][2:])
                mts.append(Mt(id=mt_id, name=words[1]))
            else:
                unparsed += option + '\n'
        return Materials(mats=mats, mts=mts, unparsed=unparsed)

    @staticmethod
    def __parse_single_mat(mat_card):
        opts = mat_card.split()
        mat = Material(mat_id=int(opts[0][1:]))
        i = 1
        while i < len(opts)-1:
            try:
                value = float(opts[i + 1])
                mat.update_nuclide(Nuclide(opts[i], opts[i + 1]))
                i += 2
            except ValueError:
                print(" Warning: mat " + str(opts[0][1:]) + " was not fully parsed!")
                break
        return mat

    @staticmethod
    def __parse_surface_macrobody(content):
        surfaces = []
        macrobodys = []
        unparsed = ''
        for option in content:
            surf_boundary = None
            surf_unparsed = ''
            reflect_match = re.match(r'\*', option)
            if reflect_match:
                surf_boundary = 1  # reflective boundary
                option = option[1:]
            surf_match = re.match(r'([0-9]+) ', option)
            if surf_match:
                words = option.split()
                surf_pair = None
                surf_tr = None  # transformation id for the surf
                surf_id = int(words[0])
                if re.match(r'[0-9\+\-]+', words[1]):
                    var = int(words[1])
                    if var < 0:
                        surf_boundary = 3  # 周期边界条件
                        surf_pair = abs(var)
                    if var > 0:
                        surf_tr = abs(var)
                    surf_type = words[2].upper()
                    other_vars = ' '.join(words[2:])
                else:
                    surf_type = words[1].upper()
                    other_vars = ' '.join(words[1:])
                if surf_type in Surface.surf_type_para:  # the surface case
                    [surf, unpar] = PlainParser._parse_option(other_vars, Surface.surf_type_para)
                    if unpar is not '':
                        surf_unparsed += ' Warning: No parsed card " ' + str(unpar) + ' " in surf ' + str(surf_id)
                    surface = Surface(number=surf_id, stype=surf_type, parameters=surf[surf_type],
                                      boundary=surf_boundary, pair=surf_pair, tr=Transformation(num=surf_tr),
                                      unparsed=surf_unparsed)
                    surfaces.append(surface)
                elif surf_type in MacroBody.body_type_para:  # the macrobody case
                    [body, unpar] = PlainParser._parse_option(other_vars, MacroBody.body_type_para)
                    if unpar is not '':
                        surf_unparsed += ' Warning: No parsed card " ' + str(unpar) + ' " in body ' + str(surf_id)
                    body = MacroBody(number=surf_id, type=surf_type, params=body[surf_type],
                                     tr=Transformation(num=surf_tr), unparsed=surf_unparsed)
                    macrobodys.append(body)
            else:
                unparsed += option + '\n'

        return Surfaces(surfaces=surfaces, macrobodys=macrobodys, unparsed=unparsed)

    @staticmethod
    def __parse_tr(content):
        trs = []
        unparsed = ''
        for option in content:
            tr_angle = None
            angle_match = re.match(r'\*', option)
            if angle_match:
                tr_angle = 1  # trn card defined by angles
                option = option[1:]
            tr_match = re.match(r'tr', option, re.I)
            if tr_match:
                words = option.split()
                num = int(words[0][2:])
                [tr_paras, index] = PlainParser._parse_list(words, float)
                tr = Transformation(paras=tr_paras, num=num, angle=tr_angle)
                tr.process()
                trs.append(tr)
            else:
                unparsed += option + '\n'
        return trs

    @staticmethod
    def __parse_geometry(content):
        cells = []
        universes = []
        geo_unparsed = ''
        for cell in content:
            cell_len = len(cell.split())
            cell_id = int(cell.split()[0])
            unparsed = ''

            if cell.split()[1].upper() != 'LIKE':  # like 'j m d geom params'
                index = 0
                # 解析几何中的材料信息
                mat_id = int(cell.split()[1])
                mat_density = None
                if mat_id == 0:
                    index = 2
                else:
                    index = 3
                    mat_density = float(cell.split()[2])

                # 解析几何中的面信息
                cell_geom = ''
                geom_no_end = True
                while geom_no_end:
                    options = cell.split()[index:]
                    for param in Cell.card_option_types:
                        if re.match(param, options[0], re.I):
                            geom_no_end = False
                            break
                    if not geom_no_end:
                        break
                    if options[0] in [':', '(', ')']:
                        cell_geom = cell_geom[0:len(cell_geom) - 1]
                        cell_geom += options[0] + ' '
                    else:
                        cell_geom += '(' + options[0] + ')' + '&'
                    index += 1
                    if index >= cell_len:
                        geom_no_end = False
                        break

                cell_geom = cell_geom[0:len(cell_geom) - 1].replace(' ', '')
                cell_geom = cell_geom.replace(')(', ')&(')
                cell_geom = cell_geom.replace('()&', '')

                # 解析几何中的其他选项
                cell_dict = {}
                cell_card_options = copy.deepcopy(Cell.card_option_types)
                if index != cell_len:
                    unparsed_items = ' '.join(cell.split()[index:])
                    while unparsed_items is not '':
                        if 'LAT' in cell_dict and 'FILL' in cell_card_options:
                            cell_card_options['FILL'] = ['list', int, -1]
                        [cell_dict_item, unparsed_items_new] = PlainParser._parse_option(unparsed_items,
                                                                                         cell_card_options)
                        if cell_dict_item:
                            for key in cell_dict_item:
                                cell_dict[key] = cell_dict_item[key]
                        if unparsed_items_new is '':
                            unparsed_items = unparsed_items_new
                        elif unparsed_items.split()[0] == unparsed_items_new.split()[0]:
                            unparsed += ' ' + unparsed_items.split()[0]
                            unparsed_items = ' '.join(unparsed_items.split()[1:])
                        else:
                            unparsed_items = unparsed_items_new

                parsed_cell = Cell(number=cell_id, material=mat_id, density=mat_density, bounds=cell_geom,
                                   unparsed=unparsed)
                parsed_cell.add_options(cell_dict)
                if parsed_cell.trcl is not None:
                    parsed_cell.trcl.process()
                cells.append(parsed_cell)

            else:  # like 'j LIKE n BUT list'
                likeid = int(cell.split()[2])
                unparsed_items = ' '.join(cell.split()[4:])
                while unparsed_items is not '':
                    if 'LAT' in cell_dict and 'FILL' in cell_card_options:
                        cell_card_options.pop('FILL')
                    [cell_dict_item, unparsed_items_new] = PlainParser._parse_option(unparsed_items,
                                                                                     cell_card_options)
                    if cell_dict_item:
                        for key in cell_dict_item:
                            cell_dict[key] = cell_dict_item[key]
                    if unparsed_items_new is '':
                        unparsed_items = unparsed_items_new
                    elif unparsed_items.split()[0] == unparsed_items_new.split()[0]:
                        unparsed += ' ' + unparsed_items.split()[0]
                        unparsed_items = ' '.join(unparsed_items.split()[1:])
                    else:
                        unparsed_items = unparsed_items_new
                # geo_unparsed += cell + '\n'

                parsed_cell = Cell(number=cell_id, likeid=likeid, unparsed=unparsed)
                parsed_cell.add_options(cell_dict)
                if parsed_cell.trcl is not None:
                    parsed_cell.trcl.process()
                for cell in cells:
                    if cell.number == likeid:
                        if parsed_cell.material is None:
                            parsed_cell.material = cell.material
                        if parsed_cell.bounds is '':
                            parsed_cell.bounds = cell.bounds
                        if parsed_cell.density is None:
                            parsed_cell.density = cell.density
                        if parsed_cell.universe is None:
                            parsed_cell.universe = cell.universe
                        if parsed_cell.lat is None:
                            parsed_cell.lat = cell.lat
                        if parsed_cell.fill is None:
                            parsed_cell.fill = cell.fill
                        break
                cells.append(parsed_cell)
        return [cells, geo_unparsed]

    @staticmethod
    def __parse_sdef(content):
        source_list = []
        distribution_list = []
        unparsed = ''
        for option in content:
            new_option = []
            for value in option.split(' '):
                if re.match(r'D\d+', value, flags=re.I):
                    value_list = list(value)
                    for final_value in value_list:
                        new_option.append(final_value)
                else:
                    new_option.append(value)
            option = ' '.join([str(x) for x in new_option])
            if option.split()[0].upper() == 'SDEF':
                source = Source()
                opt_lst = PlainParser._parse_multioption(option.split(' ', 1)[1], Source.card_option_types)
                source.add_options(opt_lst)
                source_list.append(source)
            elif re.match(r'SP', option, flags=re.I):
                distribution = Distribution()
                sp_list = PlainParser._parse_multioption(option, Distribution.card_option_types)
                distribution.add_options(sp_list)
                distribution_list.append(distribution)
            elif re.match(r'SI', option, flags=re.I):
                distribution = Distribution()
                si_list = PlainParser._parse_multioption(option, Distribution.card_option_types)
                distribution.add_options(si_list)
                distribution_list.append(distribution)
            elif re.match(r'SB', option, flags=re.I):
                distribution = Distribution()
                sb_list = PlainParser._parse_multioption(option, Distribution.card_option_types)
                distribution.add_options(sb_list)
                distribution_list.append(distribution)
            elif re.match(r'DS', option, flags=re.I):
                distribution = Distribution()
                ds_list = PlainParser._parse_multioption(option, Distribution.card_option_types)
                distribution.add_options(ds_list)
                distribution_list.append(distribution)
            else:
                unparsed += option
        return ExternalSource(source=source_list, distributions=distribution_list, unparsed=unparsed)

    @staticmethod
    def __parse_critical(content):
        kcode_list = []
        ksrc_list = []
        unparsed = ''
        for option in content:
            if option.split()[0].upper() == 'KCODE':
                kcode_list = PlainParser._parse_multioption(option, Criticality.card_option_types)
            elif option.split()[0].upper() == 'KSRC':
                ksrc_list = PlainParser._parse_multioption(option, Criticality.card_option_types)
            else:
                unparsed += option
        return Criticality(kcode=kcode_list, ksrc=ksrc_list, unparsed=unparsed)

    @staticmethod
    def _parse_option(content, cards):
        content = content.replace(' = ', ' ')
        options = content.split()
        options_dict = {}
        unparsed = []
        index = 0
        options[0] = options[0].upper()

        matched_key = None
        for key in cards.keys():
            if re.match(key, options[0], re.I) and key[-1] == options[0][-1]:
                matched_key = key
                break
        if matched_key is None:
            unparsed = ' '.join(options[index:])
            return [options_dict, unparsed]
        else:
            dtype = cards[matched_key]
            if dtype[0] == 'list':
                [opt_val, index] = PlainParser._parse_list(options, dtype[1])
            else:
                [opt_val, index] = PlainParser._parse_val(options, dtype[0])

            options_dict[matched_key] = opt_val
            # 移除已经解析过的部分
            unparsed = ' '.join(options[index:])
            return [options_dict, unparsed]

    @staticmethod
    def _parse_multioption(content, cards):
        options = content.replace(' = ', ' ').split()
        options_dict = {}
        while options:
            options[0] = options[0].upper()
            if options[0] in cards:
                dtype = cards[options[0]]
                if dtype[0] == 'list':
                    [opt_val, index] = PlainParser._parse_multilist(options, dtype[1])
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
        if options[1] == "(":
            options.remove("(")
            options.remove(")")
        index = 1
        opt_list = []
        while index < len(options):
            try:
                # to parse the fill card
                if options[index] == ":":
                    opt_list.append(options[index])
                    index += 1
                    continue
                value = func(options[index])
                index += 1
                opt_list.append(value)
            except ValueError:
                break
        # todo: check
        return [opt_list, index]

    @staticmethod
    def _parse_multilist(options, func=str):
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
                if options[0].upper() != 'SP':
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
            elif re.fullmatch(r'FCEL', options[index].upper()):
                option = options[index]
                opt_list.append(option)
                opt_list.append('=')
                index += 1
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
        else:
            return [func(options[1]), 2]
