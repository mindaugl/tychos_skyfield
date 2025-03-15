"""
Library containing Tychos object interface with Skyfield.
"""

from skyfield.vectorlib import VectorFunction
from numpy import zeros
from tychos_skyfield.baselib import TychosSystem

class ReferencePlanet():
    """
    Class for reference objects definition - used for getting both Tychos and Skyfield
    reference objects.
    ref_name: string, Tychos reference name
    skyfield_objs: skyfield.jpllib.SpiceKernel, container of Skyfield objects
    skyfield_name: Optional[string] = None, Skyfield reference name to be used directly

    Attributes
    ----------
    name
    skyfield_objs
    skyfield_name
    skyfield_obj

    """
    def __init__(self, name, skyfield_objs, skyfield_name = None):
        self.name = name.lower()
        self.skyfield_objs = skyfield_objs
        self.skyfield_name = skyfield_name
        self.skyfield_obj = self._get_skyfield_obj()

    def _get_skyfield_obj(self):
        """
        Get Skyfield object to be used as reference in TychosSkyfield
        :return: skyfield.vectorlib.VectorSum
            Skyfield object associated with the ref_name (or skyfield_name)
        """
        if self.skyfield_name is None:
            skyfield_names_dict = self.skyfield_objs.names()
            skyfield_names = list(skyfield_names_dict.values())
            skyfield_names = [x.lower() for xs in skyfield_names for x in xs]
            if self.name in skyfield_names:
                self.skyfield_name = self.name
            elif f"{self.name}_barycenter" in skyfield_names:
                self.skyfield_name = f"{self.name}_barycenter"
            else:
                raise AttributeError(f"There does not seem to be skyfield object associated "
                                     f"with the tychos object {self.name}. "
                                     f"All available skyfield objects: {skyfield_names}")
        return self.skyfield_objs[self.skyfield_name]


class TychosSkyfield(VectorFunction, TychosSystem):
    """
    Class that can be used to create tychos object in Skyfield that can be as Skyfield
    native object.
    To initialize the class, need to provide:
        - name: string for the object name in Tychos
        - ref_name: string name for the tychos object which to use as reference
        - ref_obj: ReferencePlanet object for reference name in Tychos and associated
        object in Skyfield

    Attributes
    ----------
    center
    target
    name
    ref_obj

    Methods
    -------
    at
    native_object

    NOTE: Some Skyfield routines make multiple calls for _.at() for different times.
    At the end the state will be at last called time (can be checked via self.julian_day)
    """

    def __init__(self, name, ref_obj, center=0):
        self.center = center
        self.name = name.lower()
        self.ref_obj = ref_obj
        self.target = self._get_target()
        super().__init__()

    def _get_target(self):
        """
        Get the target code to be consistent with Skyfield Ephemeris
        :return: Int or None
        """
        skyfield_names_dict = self.ref_obj.skyfield_objs.names()
        for k, v in skyfield_names_dict.items():
            if self.name == v[-1].lower():
                return k
        for k, v in skyfield_names_dict.items():
            if self.name in v[-1].lower():
                return k
        return None

    def _at(self, t):
        """
        Evaluate relative position to the tychos ref object and add skyfield ref object position.
        NOTE: Velocity of object is set to 0 (which is not correct, but does not
        affect most calculations)
        :param t: skyfield.timelib.Time
            The Time to which move the Tychos system
        :return: tuple[ndarray, ndarray, None, None]
            first element corresponds to relative position with skyfield ref object position added,
            the rest of elements to comply with the skyfield infrastructure
        """

        self.move_system(t.tt)
        obj = self[self.name]
        p = obj.location_transformed(self[self.ref_obj.name], None)
        p += self.ref_obj.skyfield_obj.at(t).position.au
        v = zeros(p.shape)
        return p, v, None, None

    def at(self, t):
        """
        Override at() method to correctly get relative position
        :param t: skyfield.timelib.Time
            The Time to which move the Tychos system
        :return: skyfield.positionlib.ICRF
            Updated position object
        """

        pos = super().at(t)
        pos.position.au -= self.ref_obj.skyfield_obj.at(t).position.au
        return pos

    def native_object(self, object_name=None):
        """
        Returns native PlanetObj associated with this object (name).
        Mostly useful for debugging.
        :param object_name: Optional[string] = None
            the name of the native object to return, default corresponds to the name of this object
        :return: tychosbaselib.PlanetObj
            native object
        """

        if object_name is None:
            object_name = self.name
        return super().__getitem__(object_name)


class Ephemeris():
    """
    Ephemeris object for Tychos system.
    For ref object name (typically 'earth'), returns Skyfield object (and not Tychos object).
    Ephemeris attribute of ref object (object.ephemeris) points to Skyfield ephemeris.
    'ref_obj' is of type ReferencePlanet and 'item' is of type string.
    Usage:
        eph_t = EphemerisTychos(ref_obj)
        eph_t['jupiter'].at(time)
    """

    def __init__(self, ref_obj):
        self.ref = ref_obj
    def __getitem__(self, item):
        name = item.lower()
        if name == self.ref.name:
            return self.ref.skyfield_obj
        return TychosSkyfield(name, self.ref)
