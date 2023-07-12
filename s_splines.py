import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Union, Optional, Dict, Set, Any, Type
import numpy as np

class ProfileParams:
    """
    Specigies the parameters to which the velocity profile must adhere
    """
    def __init__(self, dist : float, v_init : float, v_final : float, a_max : float, j_max : float, delta_t : float):
        self.dist : float = dist
        self.v_init : float = v_init
        self.v_final : float = v_final
        self.a_max : float = a_max
        self.j_max : float = j_max
        self.delta_t : float = delta_t


class Time:
    """
    Specifies the time of each segement
    """
    def __init__(self, t_1 : float, t_4 : float, t_cv : float):
        self.t_1 : float = t_1
        self.t_4 : float= t_4
        self.t_cv : float = t_cv

    @classmethod
    def from_profile_params(cls: Type['Time'], profile_params : ProfileParams ) -> 'Time':
        """
        Determines which type of profile needs to calculated 
        and returns the corresponding time parameters
        """
        return cls.from_profile_params_type_one(profile_params)

    # Initialize from ProfileParams
    @classmethod
    def from_profile_params_type_one(cls: Type['Time'], profile_params: 'ProfileParams') -> 'Time':
        """
        Generates the time parameters according to 
        profile type number one
        """
        t_1: float = profile_params.a_max / (profile_params.j_max * (4.0/np.pi + 1))
        t_4: float = (profile_params.v_final - profile_params.v_init) / profile_params.a_max - 3 * t_1
        d_v_max: float = profile_params.a_max * (9.0 * (t_1)**2 + 4.5 * t_1 * t_4 + 0.5 * (t_4)**2) + profile_params.v_init * (6 * t_1 + t_4)
        t_cv: float = (profile_params.dist - d_v_max) / profile_params.v_final
        return cls(t_1, t_4, t_cv)
    
    # Create them for the other Types as well


class Profile:
    """
    Calculates the complete velocity profile
    """
    def __init__(self, profile_params : ProfileParams):
        time : Time = Time.from_profile_params(profile_params) 
        self.k_1 : np.ndarray =  self.profile_time_segement_0(time, profile_params)


    def profile_time_segement_0(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        """
        Calculates the velocities within segment one
        """
        segment_time : float = time.t_1
        j_max : float = profile_params.j_max
        delta_t : float = profile_params.delta_t
        profile : np.ndarray = np.array([])

        for t in range(int(segment_time/delta_t)):
            current_time : float = i / delta_t

            jerk = np.sin( (np.pi*(t - 0.0)) /(2.0*segment_time) )


        return profile

    def get_velocity_at_time(current_time : float) -> float:
        """
        Retrieves the velocity at a specified time
        """
        return current_time