import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Union, Optional, Dict, Set, Any, Type

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
        self.k_2 : np.ndarray = self.profile_time_segement_1(time, profile_params)

    ########### Segment 0
    def profile_time_segement_0(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        """
        Calculates the velocities within segment one
        """
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_0_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_0_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        """
        Calculates the jerk of the zeroth segment
        """
        passed_time = 0.0
        return profile_params.j_max * np.sin((np.pi*(current_time - passed_time))/(2.0 * time.t_1 ))
    
    def _segment_0_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        """
        Calculates the acceleration of the zeroth segment
        """
        constant_a : float = 0.63662 * time.t_1
        return -0.63662 * time.t_1 * profile_params.j_max * np.cos((1.5708 * current_time) / time.t_1) + constant_a
    
    def _segment_0_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        return time.t_1 * (0.63662 * current_time - 0.405284 * time.t_1 * profile_params.j_max *  np.sin((1.5708 * current_time) / time.t_1)) + profile_params.v_init


    ########### Segment 1
    def profile_time_segement_1(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        """
        Calculates the velocities within segment one
        """
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_1_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_1_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        """
        Calculates the jerk of the zeroth segment
        """
        return profile_params.j_max 
    
    def _segment_1_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        """
        Calculates the acceleration of the zeroth segment
        """
        constant_a : float = self._segment_0_acc(time, profile_params, time.t_1)
        return profile_params.j_max * current_time + constant_a
    
    def _segment_1_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        constant_a : float = self._segment_0_acc(time, profile_params, time.t_1)
        constant_v : float = self._segment_0_vel(time,profile_params, time.t_1) 
        return  0.5* profile_params.j_max * current_time**2 + constant_a*current_time + constant_v

    ########### Segment 1
    def profile_time_segement_1(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        """
        Calculates the velocities within segment one
        """
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_1_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_1_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        """
        Calculates the jerk of the zeroth segment
        """
        return profile_params.j_max 
    
    def _segment_1_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        """
        Calculates the acceleration of the zeroth segment
        """
        constant_a : float = self._segment_0_acc(time, profile_params, time.t_1)
        return profile_params.j_max * current_time + constant_a
    
    def _segment_1_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        constant_a : float = self._segment_0_acc(time, profile_params, time.t_1)
        constant_v : float = self._segment_0_vel(time,profile_params, time.t_1) 
        return  0.5* profile_params.j_max * current_time**2 + constant_a*current_time + constant_v
    def get_velocity_at_time(self,current_time : float) -> float:
        """
        Retrieves the velocity at a specified time
        """
        return current_time
    

def main():
    profile_params : ProfileParams = ProfileParams(100.0,1.0,1.0,1,0.2,0.01)
    profile : Profile = Profile(profile_params)

    total_prof : np.ndarray = np.concatenate([profile.k_1,profile.k_2])

    plt.plot(range(len(total_prof)), total_prof)
    plt.show()

if __name__ == "__main__":
    main()