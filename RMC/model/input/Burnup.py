# -*- coding: utf-8 -*-
# author: Xiaoyu Guo

from RMC.model.input.base import YMLModelObject as BaseModel


class Burnup(BaseModel):
    """ RMC Burnup card
    reference: https://thu-real.pages.reallab.org.cn/RMC/docs/usersguide/%E7%87%83%E8%80%97%E8%AE%A1%E7%AE%97.html

    Attributes
    ----------
    burn_cell: list of int
        burnable cells
    activation_cell: list of int
        burnable cells for activation calculation
    time_step: list of float
        time steps for burnup calculation, units: day
    burnup_step: list of float
        burnup steps for burnup calculation, units: MWd/tHM
    power: list of float
        total power for burnup calculation, units: MW
    power_den: list of float
        power density for burnup calculation, units: W/gHM
    source_intensity: list of float
        neutron source intensity for activation calculation
    step_number: int
        total burnup step numbers
    current_number: int
        current burnup step for one-step burnup calculaiton
    single_step: bool
        the option for one-step burnup calculation
    substep: int
        substep for point burnup, used to improve numerical stability
    inherent: float, float
        filter for nuclides from point burnup to transport calculation.
        One is absorption fraction(default: 0.9999), and the other is density fraction(default:0.999)
    acelib: list of string
        nuclide temperature for ace lib
    strategy: dict
        burnup strategy, including bos, pcn, pcrr, and hspc.Meanwhile, the substep is also supported
        by turning the substep option on (substep=1)
    solver: int
        solver of point burnup equation
    parallel: bool
        parallel burnup calculaiton option
    succession: dict
        burnup succession calculation options, including singlestep, readnuc, writenuc, fmf, cumulativetime,
        cumulativeburnup, cellcmltvburnup
    depletemode: int
        burnup mode, 1: const flux, 2: const power
    xeequilibrium: dict
        xenon and iodine equilibrium options
    outputcell: list of string
        the expansions of cell whose nuclide density is required to be printed
    varymat: dict
        vary material options
    varysurf: dict
        vary surface options
    merge: dict
        random medium merge options for burnup calculation
    naecf: float
        neutron average energy
    impnuc: list of int
        the mandartory nuclides which are added to material
    fixednuc: int
        the max number of nuclides added to materials
    reacitons: dict
        nuclide reactions controller
    decaysource: dict
        gamma decaysource controller
    second_source_step:
        the burnup step where the second source calculation is executed
    unparsed: string

    """
    card_option_types = {
        'BURNCELL': ['list', int, -1],
        'ACTIVATIONCELL': ['list', int, -1],
        'TIMESTEP': ['list', float, -1],
        'BURNUPSTEP': ['list', float, -1],
        'POWER': ['list', float, -1],
        'POWERDEN': ['list', float, -1],
        'SOURCEINTENSITY': ['list', float, -1],
        'SUBSTEP': [int],
        'INHERENT': [float],
        'ACELIB': ['list', str, -1],
        'STRATEGY': {
            'TYPE': [int],
            'SUBSTEP': [bool],
            'SUBSTEPNUMBER': [int]
        },
        'SOLVER': [int],
        'PARALLEL': [int],
        'DEPLETEMODE': [int],
        'OUTPUTCELL': ['list', str, -1],
        'XEEQUILIBRIUM': {
            'SKIP': [int],
            'BATCHSIZE': [int],
            'POWERMODE': [int],
            'MINRATIO': [float]
        },
        'SUCCESSION': {
            'POINTBURNUP': [bool],
            'SINGLESTEP': [bool],
            'FMF': [float],
            'CUMULATIVETIME': [float],
            'CUMULATIVEBURNUP': [float],
            'CELLCMLTVBURNUP': [bool]
        },
        'VARYMAT': {
            'STEP': [int],
            'MAT': ['list', int, -1],
            'NEWMAT': ['list', int, -1]
        },
        'VARYSURF': {
            'STEP': [int],
            'SURF': ['list', int, -1],
            'NEWSURF': ['list', int, -1]
        },
        'MERGE': {
            'LEVEL': [int],
            'UNIV': ['list', int, -1]
        },
        'NAECF': [float],
        'IMPNUC': ['list', int, -1],
        'FIXEDNUC': [int],
        'REACTIONS': {
            'NUCLIDE': [int],
            'MTS': ['list', int, -1]
        },
        'DECAYSOURCE': {
            'TYPE': [int],
            'CELL': ['list', int, -1],
            'ENERGY': ['list', float, -1],
            'SECONDSOURCESTARTSTEP': ['list', int, -1]
        }

    }

    def __init__(self, burn_cell=None, activation_cell=None, time_step=None, burnup_step=None, power=None,
                 power_den=None, source_intensity=None, substep=None, inherent=None, acelib=None, strategy=None,
                 solver=None, parallel=None, succession=None, deplete_mode=None, xeequilibrium=None, outputcell=None,
                 varymat=None, varysurf=None, merge=None, naecf=None, impnuc=None, fixednuc=None, reactions=None,
                 decay_source=None,
                 unparsed=''):

        self._burn_cell = burn_cell
        self._activation_cell = activation_cell
        self._time_step = time_step
        self._burnup_step = burnup_step
        self._power = power
        self._power_den = power_den
        self._source_intensity = source_intensity
        self._succession = succession
        self._substep = substep
        self._inherent = inherent
        self._acelib = acelib
        self._strategy = strategy
        self._solver = solver
        self._parallel = parallel
        self._burnup_xeequilibrium = xeequilibrium
        self._outputcell = outputcell
        self._deplete_mode = deplete_mode
        self._merge = merge
        self._varymat = varymat
        self._varysurf = varysurf
        self._naecf = naecf
        self._impnuc = impnuc
        self._fixednuc = fixednuc
        self._reactions = reactions
        self._decay_source = decay_source
        self._second_source_step = []

        self._single_step = 0
        if self._succession is not None:
            if 'SINGLESTEP' in self._succession:
                self._single_step = self._succession['SINGLESTEP']
        # 未解析的其他选项的字符串
        self._unparsed = unparsed

        self._current_step = 0

    @property
    def burn_cell(self):
        return self._burn_cell

    @burn_cell.setter
    def burn_cell(self, burn_cell):
        self._burn_cell = burn_cell

    @property
    def activation_cell(self):
        return self._activation_cell

    @activation_cell.setter
    def activation_cell(self, cell):
        self._activation_cell = cell

    @property
    def step_number(self):
        time_size = len(self._time_step) if self._time_step else len(self._burnup_step)
        return time_size

    @property
    def single_step(self):
        return self._single_step

    @single_step.setter
    def single_step(self, single_step):
        self._single_step = single_step

    @property
    def power(self):
        return self._power

    @power.setter
    def power(self, p):
        self._power = p

    @property
    def power_den(self):
        return self._power_den

    @power_den.setter
    def power_den(self, pd):
        self._power_den = pd

    @property
    def source_intensity(self):
        return self._source_intensity

    @source_intensity.setter
    def source_intensity(self, i):
        self._source_intensity = i

    @property
    def succession(self):
        return self._succession

    @succession.setter
    def succession(self, succession):
        self._succession = succession

    def set_succession_entry(self, key, value):
        if key not in self._succession:
            raise KeyError("The key {} does not exist in burnup succession".format(key))
        self._succession[key] = value

    def del_succession_entry(self, key):
        if key not in self._succession:
            return
        del self._succession[key]

    @property
    def time_step(self):
        if self._time_step is None:
            self._time_step = []
        return self._time_step

    @time_step.setter
    def time_step(self, t):
        self._time_step = t

    @property
    def burnup_step(self):
        if self._burnup_step is None:
            self._burnup_step = []
        return self._burnup_step

    @burnup_step.setter
    def burnup_step(self, bs):
        self._burnup_step = bs

    @property
    def substep(self):
        return self._substep

    @substep.setter
    def substep(self, sub_step):
        self._substep = sub_step

    @property
    def inherent(self):
        return self._inherent

    @inherent.setter
    def inherent(self, inherent):
        self._inherent = inherent

    @property
    def acelib(self):
        return self._acelib

    @acelib.setter
    def acelib(self, ace_lib):
        self._acelib = ace_lib

    @property
    def strategy(self):
        return self._strategy

    @strategy.setter
    def strategy(self, strategy):
        self._strategy = strategy

    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, solver):
        self._solver = solver

    @property
    def parallel(self):
        return self._parallel

    @parallel.setter
    def parallel(self, parallel):
        self._parallel = parallel

    @property
    def deplete_mode(self):
        return self._deplete_mode

    @deplete_mode.setter
    def deplete_mode(self, m):
        self._deplete_mode = m

    @property
    def burnup_xeequilibrium(self):
        return self._burnup_xeequilibrium

    @burnup_xeequilibrium.setter
    def burnup_xeequilibrium(self, xeequilibrium):
        self._burnup_xeequilibrium = xeequilibrium

    @property
    def outputcell(self):
        return self._outputcell

    @outputcell.setter
    def outputcell(self, outputcell):
        self._outputcell = outputcell

    @property
    def merge(self):
        return self._merge

    @merge.setter
    def merge(self, merge):
        self._merge = merge

    @property
    def varymat(self):
        return self._varymat

    @varymat.setter
    def varymat(self, vary_mat):
        self._varymat = vary_mat

    @property
    def varysurf(self):
        return self._varysurf

    @varysurf.setter
    def varysurf(self, vary_surf):
        self._varysurf = vary_surf

    @property
    def naecf(self):
        return self._naecf

    @naecf.setter
    def naecf(self, e):
        self._naecf = e

    @property
    def impnuc(self):
        return self._impnuc

    @impnuc.setter
    def impnuc(self, nuc):
        self._impnuc = nuc

    @property
    def fixednuc(self):
        return self._fixednuc

    @fixednuc.setter
    def fixednuc(self, num):
        self._fixednuc = num

    @property
    def reactions(self):
        return self._reactions

    @reactions.setter
    def reactions(self, r):
        self._reactions = r

    @property
    def decay_source(self):
        return self._decay_source

    @decay_source.setter
    def decay_source(self, s):
        self._decay_source = s
    
    @property
    def second_source_step(self):
        if self._decay_source and 'SECONDSOURCESTARTSTEP' in self._decay_source:
            self._second_source_step = self._decay_source['SECONDSOURCESTARTSTEP']
        return self._second_source_step

    def check(self):

        if self._time_step is [] and self._burnup_step is []:
            raise TypeError("one of TIMESTEP or BURNUPSTEP need to be provided")
        if self._time_step and self._burnup_step:
            raise TypeError("TIMESTEP and BURNUPSTEP are repeatly defined in BURNUP block")
        if not ((self._power and not self._power_den and not self._source_intensity) or
                (not self._power and self._power_den and not self._source_intensity) or
                (not self._power and not self._power_den and self._source_intensity)):
            raise TypeError(
                "One of POWER, POWERDEN, and SOURCERATE need to be provided and they can't be defined simultaneously!\n")

        time_size = len(self._time_step) if self._time_step is not None else len(self._burnup_step)
        power_size = max(len(self._power), len(self._power_den), len(self._source_intensity))

        if time_size != power_size:
            raise ValueError("The number of POWER/POWERDEN is not equal to that of TIMESTEP/BURNUPSTEP")

    # 根据所有的燃耗栅元，检查燃耗栅元是否位于当前模型内，如果不存在，则删除该cell
    # 对于燃耗接续计算，则需要用户在断点处补全所有的燃耗栅元
    def check_burn_cell(self, universe, init_burn_cell):
        for cell_id in init_burn_cell:
            cell_num = universe.count_cell(cell_id=cell_id)
            # 如果模型栅元不在模型内，且燃耗栅元中有该cell_id，则删除该cell_id
            if cell_num == 0 and cell_id in self._burn_cell:
                self._burn_cell.remove(cell_id)
            # cell_id在模型内且cell_id不在burncell中，添加该cell_id
            elif cell_num > 0 and cell_id not in self._burn_cell:
                self._burn_cell.append(cell_id)
            else:
                continue

    def __str__(self):
        card = 'BURNUP\nBurnCell '
        card += ' '.join([str(x) for x in self._burn_cell]) + '\n'
        if self._activation_cell:
            card += 'ActivationCell '
            card += ' '.join([str(x) for x in self._activation_cell]) + '\n'
        if self._time_step:
            card += 'TimeStep '
            card += ' '.join([str(x) for x in self._time_step]) + '\n'
        if self._burnup_step:
            card += 'BurnupStep '
            card += ' '.join([str(x) for x in self._burnup_step]) + '\n'
        if self._power:
            card += 'Power '
            card += ' '.join([str(x) for x in self._power]) + '\n'
        if self._power_den:
            card += 'PowerDen '
            card += ' '.join([str(x) for x in self._power_den]) + '\n'
        if self._source_intensity:
            card += 'SourceIntensity '
            card += ' '.join([str(x) for x in self._source_intensity]) + '\n'
        if self._substep is not None:
            card += 'Substep ' + str(self._substep) + '\n'
        if self._inherent is not None:
            card += 'Inherent ' + str(self._inherent) + '\n'
        if self._acelib is not None:
            card += 'Acelib ' + ' '.join([str(s) for s in self._acelib]) + '\n'
        if self._strategy is not None:
            card += 'Strategy'
            card += ' type = {}'.format(self._strategy['TYPE'])
            assert "TYPE" in self._strategy.keys(), "Type option must exist in Strategy card!\n"
            if 'SUBSTEP' in self._strategy.keys():
                card += ' substep = {:d}'.format(self._strategy['SUBSTEP'])
            if 'SUBSTEPNUMBER' in self._strategy.keys():
                card += ' substepnumber = {}'.format(self._strategy['SUBSTEPNUMBER'])
            card += '\n'
        if self._solver is not None:
            card += 'Solver ' + str(self._solver) + '\n'
        if self._parallel is not None:
            card += 'Parallel ' + str(self._parallel) + '\n'
        if self._deplete_mode is not None:
            card += 'DepleteMode ' + str(self._deplete_mode) + '\n'
        if self._outputcell is not None:
            card += 'OutputCell ' + ' '.join([str(x) for x in self._outputcell]) + '\n'
        if self._merge is not None:
            card += 'Merge'
            if 'LEVEL' in self._merge:
                card += ' level = ' + str(self._merge['LEVEL'])
            if 'UNIV' in self._merge:
                card += ' univ = ' + ' '.join([str(x) for x in self._merge['UNIV']])
            card += '\n'

        if self._succession is not None:
            card += 'Succession'
            if self._single_step:
                card += ' singlestep = 1'
            if 'POINTBURNUP' in self._succession:
                if self._succession['POINTBURNUP']:
                    card += ' pointburnup = 1'
            if 'FMF' in self._succession:
                card += ' fmf = ' + str(self._succession['FMF'])
            if 'CUMULATIVETIME' in self._succession:
                card += ' cumulativetime = ' + str(self._succession['CUMULATIVETIME'])
            if 'CUMULATIVEBURNUP' in self._succession:
                card += ' cumulativeburnup = ' + str(self._succession['CUMULATIVEBURNUP'])
            if 'CELLCMLTVBURNUP' in self._succession:
                card += ' cellcmltvburnup = 1'
            card += '\n'
        if self._varymat:
            card += Burnup._dict_str(self._varymat,
                                     'VARYMAT', 'STEP', 'MAT', 'NEWMAT')
            card += '\n'
        if self._varysurf:
            card += Burnup._dict_str(self._varysurf, 'VARYSURF', 'STEP', 'SURF', 'NEWSURF')
            card += '\n'
        if self._reactions:
            card += Burnup._dict_str(self._reactions, 'REACTIONS', 'NUCLIDE', 'MTS')
            card += '\n'
        if self._decay_source:
            card += 'DecaySource'
            if 'TYPE' in self._decay_source:
                card += ' type = ' + str(self._decay_source['TYPE'])
            if 'CELL' not in self._decay_source:
                raise TypeError("CELL need to be provided in DECAYSOURCE card!\n")
            card += ' cell = ' + ' '.join([str(x) for x in self._decay_source['CELL']])
            if 'ENERGY' in self._decay_source:
                card += ' energy = ' + ' '.join([str(x) for x in self._decay_source['ENERGY']])
            card += '\n'
        if self._naecf:
            card += 'Naecf ' + str(self._naecf) + '\n'
        if self._impnuc:
            card += 'Impnuc ' + ' '.join([str(x) for x in self._impnuc]) + '\n'
        if self._fixednuc:
            card += 'FixedNuc ' + str(self._fixednuc) + '\n'
        if self._unparsed != '':
            card += self._unparsed
        card += '\n\n'
        return card

    @staticmethod
    def _dict_str(vary_info_list, vary_key, *args):
        vary_str = vary_key + ' '
        for vary_info_dict in vary_info_list:
            for vary_title in args:
                if vary_title in vary_info_dict:
                    if Burnup.card_option_types[vary_key][vary_title][0] == 'list':
                        vary_str += ' ' + vary_title.lower() + ' ='
                        for num in vary_info_dict[vary_title]:
                            vary_str += ' ' + str(num)
                    else:
                        vary_str += ' ' + vary_title.lower() + ' = ' + \
                                    str(vary_info_dict[vary_title])
                else:
                    raise ValueError(
                        'Wrong Variable Information Title in Burnup\n')
            vary_str += ' '
        return vary_str
