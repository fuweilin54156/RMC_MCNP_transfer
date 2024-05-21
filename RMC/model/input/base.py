# -*- coding:utf-8 -*-
# author: Kaiwen Li
# date: 2019-11-16

import yaml
import abc


# todo: change all of the classes inherited from this class, from __init__ to __new__, and
#       construct a __dict__ checking process.
class YMLModelObject(yaml.YAMLObject):
    """
    The method utilized by YAML package to construct Python object:
    1. use cls.__new__(cls) to create a new object
    2. use __setstate__ to update the parameter list
       or use __dict__.update() to directly update the parameter list
    3. yield the object.
    Note that, that process will not call the __init__ method, and some unexpected parameters
    can be introduced into the parameter dict, therefore, a checking method is required.
    """

    # todo: check method should be called after construction, and the __dict__ should be checked
    #       to avoid unexpected parameters from wrong input files.
    @abc.abstractmethod
    def check(self):
        pass

    @abc.abstractmethod
    def __str__(self):
        pass

    def postprocess(self):
        pass


class Model(YMLModelObject):
    def __init__(self, model=None):
        if model is None:
            model = {
                'geometry': None,
                'surface': None,
                'macrobody': None,
                'material': None,
                'refuelling': None,
                'includematerial': None,
                'unparsed': [],
                'criticality': None,
                'burnup': None,
                'criticalitysearch': None,
                'tally': None,
                'print': None,
                'mesh': None,
                'binaryout': None,
                'externalsource': None,
                'ptrac': None,
                'physics': None,
                'fixedsource': None
            }
        self.model = model

    def check(self):
        for key in self.model.keys():
            if key is not 'unparsed' and self.model[key] is not None:
                self.model[key].check()

    def postprocess(self):
        self.check()
        if self.model['surface']:
            surf_max_id = self.model['surface'].get_max_surf_id()
        else:
            surf_max_id = 0
        if self.model['macrobody']:
            for macrobody in self.model['macrobody'].macrobodies.values():
                surfaces = macrobody.transfer(max_surf_id=surf_max_id)
                if self.model['surface']:
                    self.model['surface'].surfaces.extend(surfaces.surfaces)
                else:
                    self.model['surface'] = surfaces
                surf_max_id += len(surfaces.surfaces)
            if self.model['tally']:
                self.model['tally'].pro_tally(macrobodies=self.model['macrobody'].macrobodies)
            if self.model['binaryout']:
                self.model['binaryout'].pro_binaryout(macrobodies=self.model['macrobody'].macrobodies)
            if self.model['externalsource']:
                self.model['externalsource'].pro_source(macrobodies=self.model['macrobody'].macrobodies)
            if self.model['ptrac']:
                self.model['ptrac'].pro_ptrac(macrobodies=self.model['macrobody'].macrobodies)
            self.model['geometry'].trans_cell(macrobodies=self.model['macrobody'].macrobodies)
            del self.model['macrobody']
        post_surfs = self.model['geometry'].proc_lat5(surf_max_id)  # process the lat=5 cases
        if post_surfs:
            for surf in post_surfs:
                self.model['surface'].add_surface(surf)
        if self.model['geometry'] is not None:
            self.model['geometry'].postprocess()
        if self.model['externalsource'] is not None:
            self.model['externalsource'].postprocess()
            # process the cell expansion
            for source in self.model['externalsource'].source:
                source.postprocess(geometry=self.geometry)
            for distribution in self.model['externalsource'].distributions:
                distribution.postprocess(geometry=self.geometry)

    # todo: the sequence of output has not been defined.
    def __str__(self):
        s = ''
        for key in self.model.keys():
            if key is not 'unparsed' and self.model[key] is not None:
                test_str_=str(self.model[key])
                # print(s)
                # print("####################")
                # print(test_str_)
                s += str(self.model[key])
        for card in self.model['unparsed']:
            s += card + '\n\n'
        return s.strip('\n') + '\n'

    @property
    def geometry(self):
        return self.model['geometry']

    def __getitem__(self, item):
        if item in self.model:
            return self.model[item]
        else:
            return None

    def __setitem__(self, item, value):
        if item in self.model:
            self.model[item] = value
        else:
            raise AttributeError('Model ' + item + 'is not supported.')
