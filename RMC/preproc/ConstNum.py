import os
import RMC.preproc.Base as FB


class ConstNum:
    def __init__(self, filename):
        self.numlist = {
            'PI': 3.14159265358979323846,
            'NA': 6.02214076e+23,
            'H1': 1.00782,
            'D2': 2.0141,
            'T3': 3.01604,
            'He3': 3.01602,
            'He4': 4.0026,
            'Li6': 6.01512,
            'Li7': 7.016,
            'Be9': 9.01218,
            'B10': 10.012936,
            'B11': 11.009305,
            'C12': 12.0,
            'C13': 13.003354,
            'C14': 14.003241,
            'N14': 14.003074,
            'N15': 15.000108,
            'O16': 15.994914,
            'O17': 16.999131,
            'O18': 17.999159,
            'F19': 18.998403,
            'Ne20': 19.99244,
            'Ne21': 20.993846,
            'Ne22': 21.991385,
            'Na23': 22.989769,
            'Mg24': 23.985041,
            'Mg25': 24.985836,
            'Mg26': 25.982592,
            'Al27': 26.981538,
            'Si28': 27.976926,
            'Si29': 28.976494,
            'Si30': 29.97377,
            'P31': 30.973761,
            'S32': 31.972071,
            'S33': 32.971458,
            'S34': 33.967867,
            'S36': 35.96708,
            'Cl35': 34.968852,
            'Cl37': 36.965902,
            'Ar36': 35.967545,
            'Ar38': 37.962732,
            'Ar40': 39.962383,
            'K39': 38.963706,
            'K40': 39.963998,
            'K41': 40.961825,
            'Ca40': 39.96259,
            'Ca42': 41.958617,
            'Ca43': 42.958766,
            'Ca44': 43.955481,
            'Ca46': 45.953689,
            'Ca48': 47.952522,
            'Sc45': 44.955908,
            'Ti46': 45.952627,
            'Ti47': 46.951758,
            'Ti48': 47.947941,
            'Ti49': 48.947865,
            'Ti50': 49.944786,
            'V50': 49.947156,
            'V51': 50.943957,
            'Cr50': 49.946041,
            'Cr52': 51.940506,
            'Cr53': 52.940648,
            'Cr54': 53.938879,
            'Mn55': 54.938043,
            'Fe54': 53.939608,
            'Fe56': 55.934936,
            'Fe57': 56.935392,
            'Fe58': 57.933274,
            'Co59': 58.933194,
            'Ni58': 57.935342,
            'Ni60': 59.930785,
            'Ni61': 60.931055,
            'Ni62': 61.928345,
            'Ni64': 63.927966,
            'Cu63': 62.929597,
            'Cu65': 64.927789,
            'Zn64': 63.929142,
            'Zn66': 65.926033,
            'Zn67': 66.927127,
            'Zn68': 67.924844,
            'Zn70': 69.925319,
            'Ga69': 68.925573,
            'Ga71': 70.924702,
            'Ge70': 69.924248,
            'Ge72': 71.922075,
            'Ge73': 72.923458,
            'Ge74': 73.921177,
            'Ge76': 75.921402,
            'As75': 74.921594,
            'Se74': 73.922475,
            'Se76': 75.919213,
            'Se77': 76.919914,
            'Se78': 77.917309,
            'Se80': 79.916521,
            'Se82': 81.916699,
            'Br79': 78.918337,
            'Br81': 80.916289,
            'Kr78': 77.920364,
            'Kr80': 79.916378,
            'Kr82': 81.913482,
            'Kr83': 82.914127,
            'Kr84': 83.911497,
            'Kr86': 85.91061,
            'Rb85': 84.911789,
            'Rb87': 86.90918,
            'Sr84': 83.913419,
            'Sr86': 85.90926,
            'Sr87': 86.908877,
            'Sr88': 87.905612,
            'Y89': 88.90584,
            'Zr90': 89.904697,
            'Zr91': 90.905639,
            'Zr92': 91.905034,
            'Zr94': 93.90631,
            'Zr96': 95.908271,
            'Nb93': 92.906373,
            'Mo92': 91.906807,
            'Mo94': 93.905084,
            'Mo95': 94.905838,
            'Mo96': 95.904676,
            'Mo97': 96.906018,
            'Mo98': 97.905404,
            'Mo100': 99.907471,
            'Tc97': 96.906366,
            'Tc98': 97.907212,
            'Tc99': 98.90625,
            'Ru96': 95.90759,
            'Ru98': 97.905286,
            'Ru99': 98.905934,
            'Ru100': 99.904214,
            'Ru101': 100.9055769,
            'Ru102': 101.9043441,
            'Ru104': 103.9054275,
            'Rh103': 102.905498,
            'Pd102': 101.9056022,
            'Pd104': 103.9040305,
            'Pd105': 104.9050796,
            'Pd106': 105.9034804,
            'Pd108': 107.9038916,
            'Pd110': 109.905172,
            'Ag107': 106.9050916,
            'Ag109': 108.9047553,
            'Cd106': 105.9064599,
            'Cd108': 107.9041834,
            'Cd110': 109.903006,
            'Cd111': 110.904182,
            'Cd112': 111.902762,
            'Cd113': 112.904408,
            'Cd114': 113.903365,
            'Cd116': 115.904763,
            'In113': 112.904061,
            'In115': 114.903878,
            'Sn112': 111.904823,
            'Sn114': 113.9027827,
            'Sn115': 114.903344,
            'Sn116': 115.901742,
            'Sn117': 116.902953,
            'Sn118': 117.901606,
            'Sn119': 118.903311,
            'Sn120': 119.902201,
            'Sn122': 121.9034438,
            'Sn124': 123.9052766,
            'Sb121': 120.903812,
            'Sb123': 122.9042132,
            'Te120': 119.9040593,
            'Te122': 121.9030435,
            'Te123': 122.9042698,
            'Te124': 123.9028171,
            'Te125': 124.9044299,
            'Te126': 125.9033109,
            'Te128': 127.904461,
            'Te130': 129.906222,
            'I127': 126.9044719,
            'Xe124': 123.905892,
            'Xe126': 125.9042983,
            'Xe128': 127.903531,
            'Xe129': 128.90478,
            'Xe130': 129.903509,
            'Xe131': 130.905084,
            'Xe132': 131.904155,
            'Xe134': 133.905394,
            'Xe136': 135.907214,
            'Cs133': 132.905451,
            'Ba130': 129.9063207,
            'Ba132': 131.9050611,
            'Ba134': 133.904508,
            'Ba135': 134.905688,
            'Ba136': 135.904575,
            'Ba137': 136.905827,
            'Ba138': 137.905247,
            'La138': 137.9071149,
            'La139': 138.9063563,
            'Ce136': 135.907129,
            'Ce138': 137.905991,
            'Ce140': 139.9054431,
            'Ce142': 141.9092504,
            'Pr141': 140.9076576,
            'Nd142': 141.907729,
            'Nd143': 142.90982,
            'Nd144': 143.910093,
            'Nd145': 144.9125793,
            'Nd146': 145.9131226,
            'Nd148': 147.9168993,
            'Nd150': 149.9209022,
            'Pm145': 144.9127559,
            'Pm147': 146.915145,
            'Sm144': 143.9120065,
            'Sm147': 146.9149044,
            'Sm148': 147.9148292,
            'Sm149': 148.9171921,
            'Sm150': 149.9172829,
            'Sm152': 151.9197397,
            'Sm154': 153.9222169,
            'Eu151': 150.9198578,
            'Eu153': 152.921238,
            'Gd152': 151.9197995,
            'Gd154': 153.9208741,
            'Gd155': 154.9226305,
            'Gd156': 155.9221312,
            'Gd157': 156.9239686,
            'Gd158': 157.9241123,
            'Gd160': 159.9270624,
            'Tb159': 158.9253547,
            'Dy156': 155.9242847,
            'Dy158': 157.9244159,
            'Dy160': 159.9252046,
            'Dy161': 160.9269405,
            'Dy162': 161.9268056,
            'Dy163': 162.9287383,
            'Dy164': 163.9291819,
            'Ho165': 164.9303288,
            'Er162': 161.9287884,
            'Er164': 163.9292088,
            'Er166': 165.9302995,
            'Er167': 166.9320546,
            'Er168': 167.9323767,
            'Er170': 169.9354702,
            'Tm169': 168.9342179,
            'Yb168': 167.9338896,
            'Yb170': 169.9347664,
            'Yb171': 170.9363302,
            'Yb172': 171.9363859,
            'Yb173': 172.9382151,
            'Yb174': 173.9388664,
            'Yb176': 175.9425764,
            'Lu175': 174.9407752,
            'Lu176': 175.9426897,
            'Hf174': 173.9400461,
            'Hf176': 175.9414076,
            'Hf177': 176.9432277,
            'Hf178': 177.9437058,
            'Hf179': 178.9458232,
            'Hf180': 179.946557,
            'Ta180': 179.9474648,
            'Ta181': 180.9479958,
            'W180': 179.9467108,
            'W182': 181.948203,
            'W183': 182.950222,
            'W184': 183.95093,
            'W186': 185.9543628,
            'Re185': 184.9529545,
            'Re187': 186.9557501,
            'Os184': 183.9524885,
            'Os186': 185.953835,
            'Os187': 186.9557474,
            'Os188': 187.9558352,
            'Os189': 188.9581442,
            'Os190': 189.9584437,
            'Os192': 191.961477,
            'Ir191': 190.9605893,
            'Ir193': 192.9629216,
            'Pt190': 189.9599297,
            'Pt192': 191.9610387,
            'Pt194': 193.9626809,
            'Pt195': 194.9647917,
            'Pt196': 195.964952,
            'Pt198': 197.9678949,
            'Au197': 196.966568,
            'Hg196': 195.9658326,
            'Hg198': 197.966768,
            'Hg199': 198.96828,
            'Hg200': 199.968326,
            'Hg201': 200.970302,
            'Hg202': 201.970643,
            'Hg204': 203.973493,
            'Tl203': 202.9723446,
            'Tl205': 204.9744278,
            'Pb204': 203.973044,
            'Pb206': 205.9744657,
            'Pb207': 206.9758973,
            'Pb208': 207.9766525,
            'Bi209': 208.9803991,
            'Po209': 208.9824308,
            'Po210': 209.9828741,
            'At210': 209.9871479,
            'At211': 210.9874966,
            'Rn211': 210.9906011,
            'Rn220': 220.0113941,
            'Rn222': 222.0175782,
            'Fr223': 223.019736,
            'Ra223': 223.0185023,
            'Ra224': 224.020212,
            'Ra226': 226.0254103,
            'Ra228': 228.0310707,
            'Ac227': 227.0277523,
            'Th230': 230.0331341,
            'Th232': 232.0380558,
            'Pa231': 231.0358842,
            'U233': 233.0396355,
            'U234': 234.0409523,
            'U235': 235.0439301,
            'U236': 236.0455682,
            'U238': 238.0507884,
            'Np236': 236.04657,
            'Np237': 237.0481736,
            'Pu238': 238.0495601,
            'Pu239': 239.0521636,
            'Pu240': 240.0538138,
            'Pu241': 241.0568517,
            'Pu242': 242.0587428,
            'Pu244': 244.0642053,
            'Am241': 241.0568293,
            'Am243': 243.0613813,
            'Cm243': 243.0613893,
            'Cm244': 244.0627528,
            'Cm245': 245.0654915,
            'Cm246': 246.0672238,
            'Cm247': 247.0703541,
            'Cm248': 248.0723499,
            'Bk247': 247.0703073,
            'Bk249': 249.0749877,
            'Cf249': 249.0748539,
            'Cf250': 250.0764062,
            'Cf251': 251.0795886,
            'Cf252': 252.0816272,
            'Es252': 252.08298,
            'Fm257': 257.0951061,
            'Md258': 258.0984315,
            'Md260': 260.103653,
            'No259': 259.101031,
            'Lr262': 262.109612,
            'Rf267': 267.121796,
            'Db268': 268.125675,
            'Sg271': 271.133936,
            'Bh272': 272.138265,
            'Hs270': 270.134292,
            'Mt276': 276.151595,
            'Ds281': 281.164515,
            'Rg280': 280.165146,
            'Cn285': 285.177126,
            'Nh284': 284.178736,
            'Fl289': 289.190426,
            'Mc288': 288.192746,
            'Lv293': 293.204496,
            'Ts292': 292.207467,
            'Og294': 294.213927
        }
        self.file_name = os.path.abspath(filename)
        self.content = ''

    def repconstnum(self):
        self._read_in()
        '''
        常量替换思路：
        替换文档中出现的常量
        '''
        print(' replacing Constant Variable ...')
        for name, value in self.numlist.items():
            self.content = FB.Base.sub_string(self.content, [0, len(self.content)], name, str(value))
        print(' Constant Variable replaced')
        self._write_file()

    def _read_in(self):
        with open(self.file_name, 'r') as f:
            self.content = f.read()

    def _write_file(self):
        with open(self.file_name, 'w+') as f:
            f.write(self.content)