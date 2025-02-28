"""
Library containing Tychos object interface with Skyfield.
"""
from skyfield.vectorlib import VectorFunction
from numpy import array
from tychosbaselib import TychosSystem


class TychosSkyfield(VectorFunction, TychosSystem):
    """
    Class that can be used to create tychos object in Skyfield that can be as Skyfield
    native object.
    To initialize the class, need to provide:
        - name_tychos: string in the format of "'name'_tychos" where 'name' is the name of
        object, e.g., for jupiter, it is 'jupiter_tychos'.
        - ref_name_tychos: string for the tychos object which to use as reference, e.g.,
        'earth_tychos'
        - ref_obj_skyfield: 'skyfield.positionlib.Barycentric' object corresponding to the
        Barycentric position of reference object in Skyfield. It can be obtained, e.g.,
        for earth, by calling: skyfield.api.load('de421.bsp')['earth']

    Attributes
    ----------
    center
    name
    target
    ref_obj_skyfield

    Methods
    -------
    get_all_objects
    get_observable_objects

    NOTE: Some Skyfield routines make multiple calls for _.at() for different times.
    At the end the state will be at last called time (can be checked via self.julian_day)

    """

    def __init__(self, name_tychos, ref_name_tychos, ref_obj_skyfield, center=0):
        self.center = center
        self.target = None
        self.name = name_tychos
        self._name_native = self._tychos_native_name(name_tychos)
        self._ref_name_native = self._tychos_native_name(ref_name_tychos)
        self.ref_obj_skyfield = ref_obj_skyfield
        super().__init__()

    def __str__(self):
        return f"TychosSkyfield, all observable objects are: {self.get_observable_objects()}."

    def _tychos_native_name(self, name_string):
        """
        Take name_string which is of form 'name_tychos' and return 'name' for the pure
        Tychos library routines
        :param name_string: string
        :return: string
        """

        try:
            name, end = name_string.split("_")
            assert end.lower() == "tychos"
        except Exception as e:
            raise AttributeError(f"Unknown Tychos object {name_string}, "
                                 f"all possible objects are {self.get_all_objects()}.") from e
        return name

    def _at(self, time):
        """
        Evaluate relative position to the tychos ref object and add skyfield ref object position.
        :param time: skyfield.timelib.Time
            The Time to which move the Tychos system
        :return: tuple[ndarray, ndarray, None, None]
            first element corresponds to relative position with skyfield ref object position added,
            the rest of elements to comply with the skyfield infrastructure
        """

        self.move_system(time.tt)
        obj = self[self._name_native]
        p = obj.location_transformed(self[self._ref_name_native], None)
        p += self.ref_obj_skyfield.position.au
        v = array([0, 0, 0])
        return p, v, None, None

    def get_all_objects(self):
        """
        Get all tychos objects.
        Overrides base class function to include "_tychos" in the name
        :return: list[string]
        """

        return [x + "_tychos" for x in super().get_all_objects()]

    def get_observable_objects(self):
        """
        Get observable tychos objects.
        Overrides base class function to include "_tychos" in the name
        :return: list[string]
        """

        return [x + "_tychos" for x in super().get_observable_objects()]
