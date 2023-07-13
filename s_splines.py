import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Union, Optional, Dict, Set, Any, Type
from scipy.optimize import fsolve

class ProfileParams:
    """
    Specigies the parameters to which the velocity profile must adhere
    """
    def __init__(self, dist : float, v_init : float, v_final : float,v_max : float, a_max : float, j_max : float, delta_t : float):
        self.dist : float = dist
        self.v_init : float = v_init
        self.v_final : float = v_final
        self.v_max : float = v_max
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

    @staticmethod
    def calculate_it1a(j_max : float , v_initial : float, dist: float) -> float:
        # Define the equation to be solved
        def equation(it1a : float) -> float:
            return 9 * j_max * (4/np.pi + 1) * (it1a ** 3) + 6 * v_initial * it1a - dist

        # Use fsolve to find the root of the equation
        it1a = fsolve(equation, 0.1)  # 0.1 is the initial guess

        # Return the smallest positive root
        return it1a[it1a > 0].min()

    @staticmethod
    def calculate_it4a(a_max: float, j_max : float , v_initial : float, dist: float) -> float:
        # Define the equation to be solved
        def equation(it4a : float) -> float:
            A = 0.5*a_max
            B = (4.5 * (a_max **2 ) / (j_max * (4/np.pi + 1)) + v_initial)
            C = (9 * a_max **3)/((j_max * 4)/np.pi + j_max)**2 + (6*v_initial* a_max)/(j_max * (4/np.pi + 1)) - dist
            return A * it4a ** 2 + B * it4a + C

        # Use fsolve to find the root of the equation
        try:
            it4a = fsolve(equation, 0.1)  # 0.1 is the initial guess
            return it4a[it4a > 0].min()
        except ValueError:
            return False


    @classmethod
    def from_profile_params(cls: Type['Time'], profile_params : ProfileParams ) -> 'Time':
        return cls.from_profile_params_type_one(profile_params)

    # Initialize from ProfileParams
    @classmethod
    def from_profile_params_type_one(cls: Type['Time'], profile_params: 'ProfileParams') -> 'Time':
        """
        Generates the time parameters according to 
        profile type number one
        """
        print("from_profile_params_type_one")
        a_max : float = profile_params.a_max
        j_max : float = profile_params.j_max
        v_init : float = profile_params.v_init
        v_max : float = profile_params.v_max

        t_1: float = profile_params.a_max / (profile_params.j_max * (4.0/np.pi + 1))
        t_4: float = (profile_params.v_final - profile_params.v_init) / profile_params.a_max - 3 * t_1
        d_v_max: float = profile_params.a_max * (9.0 * (t_1)**2 + 4.5 * t_1 * t_4 + 0.5 * (t_4)**2) + profile_params.v_init * (6 * t_1 + t_4)
        v_a_max = 3 * a_max * t_1
        t_cv: float = (profile_params.dist - d_v_max) / profile_params.v_final

        print(f"{t_4}, {t_cv}")

        if (t_4 > 0 and t_cv > 0) and (v_a_max < v_max - v_init) :
                  return Time(t_1, t_4, t_cv)
                
        return cls.from_profile_params_type_two(profile_params)
    
      # Initialize from ProfileParams
    @classmethod
    def from_profile_params_type_two(cls: Type['Time'], profile_params: 'ProfileParams') -> 'Time':
        """
        Generates the time parameters according to 
        profile type number one
        """
        print("from_profile_params_type_two, fallback")

        a_max : float = profile_params.a_max
        j_max : float = profile_params.j_max
        v_init : float = profile_params.v_init
        v_max : float = profile_params.v_max


        t_1: float = profile_params.a_max / (profile_params.j_max * (4.0/np.pi + 1))
        t_4: float = Time.calculate_it4a(a_max,j_max,v_init, profile_params.dist)
        v_final : float = profile_params.a_max * (3 * t_1 + t_4) + profile_params.v_init
        
        if not Time.calculate_it4a(a_max,j_max,v_init, profile_params.dist):
            return cls.from_profile_params_type_three(profile_params)
        if v_final > v_max:
            return cls.from_profile_params_type_three(profile_params)

            


        profile_params.v_final = v_final
        d_v_max: float = profile_params.a_max * (9.0 * (t_1)**2 + 4.5 * t_1 * t_4 + 0.5 * (t_4)**2) + profile_params.v_init * (6 * t_1 + t_4)
        t_cv: float = round((profile_params.dist - d_v_max) / profile_params.v_final, 7)
        v_a_max = 3 * a_max * t_1 
        a_peak : float = j_max*(4/np.pi + 1)* t_1
        print(f"{t_4}, {t_cv}, {a_peak}")
        
        if ((t_4 > 0 and t_cv <= 0 ) and abs(a_max - a_peak) < 0.001) :
            return Time(t_1, t_4, t_cv)
                
        return cls.from_profile_params_type_three(profile_params)
    

          # Initialize from ProfileParams
    @classmethod
    def from_profile_params_type_three(cls: Type['Time'], profile_params: 'ProfileParams') -> 'Time':
        """
        Generates the time parameters according to 
        profile type number one
        """
        print("from_profile_params_type_three")
        j_max : float = profile_params.j_max
        v_init : float = profile_params.v_init
        v_max : float = profile_params.v_max
        a_max : float = profile_params.a_max

        t_1 : float = ((v_max - v_init)/(3*j_max * ((4.0/np.pi + 1) )))**0.5
        t_4: float = (profile_params.v_final - profile_params.v_init) / profile_params.a_max - 3 * t_1
        d_v_max2_acc : float = j_max* (36/np.pi + 9)* t_1 **3 + 6 * v_init * t_1
        t_cv: float = (profile_params.dist - d_v_max2_acc) / profile_params.v_max

        a_p : float = j_max * (4.0/np.pi + 1) * t_1
        print(f"{t_4}, {t_cv}")

        if t_4 <= 0 and t_cv > 0 :
            return Time(t_1, t_4, t_cv)

                
        return cls.from_profile_params_type_four(profile_params)
    
    @classmethod
    def from_profile_params_type_four(cls: Type['Time'], profile_params: 'ProfileParams') -> 'Time':
        """
        Generates the time parameters according to 
        profile type number one
        """
        t_1 : float= Time.calculate_it1a(profile_params.j_max , profile_params.v_init, profile_params.dist)
        t_4 : float= 0.0
        t_cv : float = 0.0
        v_final : float = profile_params.j_max * (12/np.pi + 3) * t_1 **2 + profile_params.v_init
        profile_params.v_final = v_final

        assert t_1 > 0, f"no rise time brother: {t_1}"

        return  Time(t_1, t_4, t_cv)
    
    
    
    # Create them for the other Types as well


class Profile:
    """
    Calculates the complete velocity profile
    """
    def __init__(self, profile_params : ProfileParams):
        time : Time = Time.from_profile_params(profile_params) 
        self.k_1_vel : np.ndarray =  self.profile_time_segement_1(time, profile_params)
        self.k_1_acc : np.ndarray =  self.profile_time_segement_1_acc(time, profile_params)
        self.k_1_jerk : np.ndarray =  self.profile_time_segement_1_jerk(time, profile_params)

        self.k_2_vel : np.ndarray = self.profile_time_segement_2(time, profile_params)
        self.k_2_acc : np.ndarray =  self.profile_time_segement_2_acc(time, profile_params)
        self.k_2_jerk : np.ndarray = self.profile_time_segement_2_jerk(time, profile_params)

        self.k_3_vel : np.ndarray = self.profile_time_segement_3(time, profile_params)
        self.k_3_acc : np.ndarray =  self.profile_time_segement_3_acc(time, profile_params)
        self.k_3_jerk : np.ndarray = self.profile_time_segement_3_jerk(time, profile_params)

        self.k_4_vel : np.ndarray = self.profile_time_segement_4(time, profile_params)
        self.k_4_acc : np.ndarray =  self.profile_time_segement_4_acc(time, profile_params)
        self.k_4_jerk : np.ndarray = self.profile_time_segement_4_jerk(time, profile_params)

        self.k_5_vel : np.ndarray = self.profile_time_segement_5(time, profile_params)
        self.k_5_acc : np.ndarray =  self.profile_time_segement_5_acc(time, profile_params)
        self.k_5_jerk : np.ndarray = self.profile_time_segement_5_jerk(time, profile_params)

        self.k_6_vel : np.ndarray = self.profile_time_segement_6(time, profile_params)
        self.k_6_acc : np.ndarray =  self.profile_time_segement_6_acc(time, profile_params)
        self.k_6_jerk : np.ndarray = self.profile_time_segement_6_jerk(time, profile_params)

        self.k_7_vel : np.ndarray = self.profile_time_segement_7(time, profile_params)
        self.k_7_acc : np.ndarray =  self.profile_time_segement_7_acc(time, profile_params)
        self.k_7_jerk : np.ndarray = self.profile_time_segement_7_jerk(time, profile_params)

        self.k_8_vel : np.ndarray = self.profile_time_segement_8(time, profile_params)
        self.k_8_acc : np.ndarray =  self.profile_time_segement_8_acc(time, profile_params)
        self.k_8_jerk : np.ndarray = self.profile_time_segement_8_jerk(time, profile_params)

        self.total_profile : np.ndarray = np.concatenate([self.k_1_vel,
                                                          self.k_2_vel,
                                                          self.k_3_vel,
                                                          self.k_4_vel,
                                                          self.k_5_vel,
                                                          self.k_6_vel,
                                                          self.k_7_vel,
                                                          self.k_8_vel])
        
        self.final_pos : float = Profile.final_distance(self.total_profile, profile_params.delta_t)

    @staticmethod
    def final_distance(profile: np.ndarray, dt : float) -> float:
        dist : float = 0.0

        for i in profile:
            dist += i * dt

        print(f"this is the final distance: {dist}")
        return dist


    ############################################################################################# Segment 1
    def profile_time_segement_1(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_1_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def profile_time_segement_1_acc(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_1_acc(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile    
    
    def profile_time_segement_1_jerk(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_1_jerk(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    
    def _segment_1_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time

        return J * np.sin( (np.pi * x) / (2 * T) )
    
    def _segment_1_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_0 = (2 * J * T) / np.pi
        
        return - (2 * J * T * np.cos( (np.pi * x) / (2 * T) )) / np.pi + a_0
    
    def _segment_1_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_0 = (2 * J * T) / np.pi
        v_0 = profile_params.v_init
        return  - (4 * J * T **2 * np.sin( (np.pi * x) / (2 * T) )) / (np.pi **2) + a_0*x + v_0


    ############################################################################################# Segment 2
    def profile_time_segement_2(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t +1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_2_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def profile_time_segement_2_acc(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_2_acc(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile  
    
    def profile_time_segement_2_jerk(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_2_jerk(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_2_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        return J
    
    def _segment_2_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_1 : float = self._segment_1_acc(time, profile_params , T)
        return J * x + a_1
    
    def _segment_2_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_1 : float = self._segment_1_acc(time, profile_params , T)
        v_1 : float = self._segment_1_vel(time, profile_params , T)
        return  0.5 * J * x **2.0 + a_1 * x + v_1 

    ############################################################################################# Segment 3
    def profile_time_segement_3(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_3_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile

    def profile_time_segement_3_acc(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_3_acc(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
      
    def profile_time_segement_3_jerk(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t +1 )):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_3_jerk(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_3_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        return J * np.sin(np.pi / 2 +  (np.pi * x) / (2 * T) )
    
    def _segment_3_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_3 : float = self._segment_2_acc(time, profile_params , T)
        return  (2 * J * T * np.sin( (np.pi * x) / (2 * T) )) / np.pi + a_3
    
    def _segment_3_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_3 : float = self._segment_2_acc(time, profile_params, time.t_1)
        v_3 : float = self._segment_2_vel(time,profile_params, time.t_1) + (4 * J * T **2)/ (np.pi **2)
        return   - (4 * J * T **2 * np.cos( (np.pi * x) / (2 * T) )) / (np.pi **2) + a_3*x + v_3
    ############################################################################################# Segment 4
    def profile_time_segement_4(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_4/profile_params.delta_t +1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_4_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile

    def profile_time_segement_4_acc(self, time : Time, profile_params : ProfileParams) -> np.ndarray:  
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_4/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_4_acc(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile


    def profile_time_segement_4_jerk(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_4/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_4_jerk(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_4_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        return 0
    
    def _segment_4_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_4 : float = self._segment_3_acc(time, profile_params , T)
        return a_4
    
    def _segment_4_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        x = 0 if x < 0 else x

        a_4 : float = self._segment_3_acc(time, profile_params , T)
        v_4 : float = self._segment_3_vel(time, profile_params , T)
        return  a_4 * x + v_4
    ############################################################################################# Segment 5
    def profile_time_segement_5(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_5_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def profile_time_segement_5_acc(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_5_acc(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def profile_time_segement_5_jerk(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_5_jerk(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_5_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T1,T4, x = profile_params.j_max, time.t_1,time.t_4, current_time
        return  - J * np.sin( (np.pi * x) / (2 * T1) )
    
    def _segment_5_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        J, T1,T4, x = profile_params.j_max, time.t_1,time.t_4, current_time
        a_5 : float = self._segment_4_acc(time,profile_params,T4) - (2 * J * T1)/np.pi 
        return (2 * J * T1 * np.cos( (np.pi * x) / (2 * T1) )) / np.pi + a_5
    
    def _segment_5_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T1,T4, x = profile_params.j_max, time.t_1,time.t_4, current_time
        a_5 : float = self._segment_4_acc(time,profile_params, T4) - (2 * J * T1)/np.pi 
        v_5 : float = self._segment_4_vel(time,profile_params, T4)
        return (4 * J * T1 **2 * np.sin( (np.pi * x) / (2 * T1) )) / (np.pi **2) + a_5*x +  v_5
    ############################################################################################# Segment 6
    def profile_time_segement_6(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_6_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def profile_time_segement_6_acc(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_6_acc(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def profile_time_segement_6_jerk(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_6_jerk(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_6_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        return -J
    
    def _segment_6_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_6 : float = self._segment_5_acc(time, profile_params , T)
        return -J * x + a_6
    
    def _segment_6_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_6 : float = self._segment_5_acc(time, profile_params , T)
        v_6 : float = self._segment_5_vel(time, profile_params , T)
        return  - 0.5 * J * x **2.0 + a_6 * x + v_6 
    ############################################################################################# Segment 7
    def profile_time_segement_7(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_7_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile

    def profile_time_segement_7_acc(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_7_acc(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
      
    def profile_time_segement_7_jerk(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_1/profile_params.delta_t +1 )):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_7_jerk(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_7_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        return - J * np.sin(np.pi / 2 +  (np.pi * x) / (2 * T) )
    
    def _segment_7_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_7 : float = self._segment_6_acc(time, profile_params , T)
        return  - (2 * J * T * np.sin( (np.pi * x) / (2 * T) )) / np.pi + a_7
    
    def _segment_7_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        a_7 : float = self._segment_6_acc(time, profile_params, time.t_1)
        v_7 : float = self._segment_6_vel(time,profile_params, time.t_1) - (4 * J * T **2)/ (np.pi **2)
        return    (4 * J * T **2 * np.cos( (np.pi * x) / (2 * T) )) / (np.pi **2) + a_7*x + v_7
    
    ############################################################################################# Segment 8
    def profile_time_segement_8(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_cv/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_8_vel(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile

    def profile_time_segement_8_acc(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_cv/profile_params.delta_t + 1)):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_8_acc(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
      
    def profile_time_segement_8_jerk(self, time : Time, profile_params : ProfileParams) -> np.ndarray:
        profile : np.ndarray = np.array([])
        for t in range(int(time.t_cv/profile_params.delta_t +1 )):
            current_time : float = t * profile_params.delta_t
            vel = self._segment_8_jerk(time, profile_params, current_time)
            profile = np.append(profile, vel)
        return profile
    
    def _segment_8_jerk(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        return 0
    
    def _segment_8_acc(self, time : Time, profile_params : ProfileParams, current_time : float, ) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        return  0
    
    def _segment_8_vel(self, time : Time,profile_params : ProfileParams, current_time : float) -> float:
        J, T, x = profile_params.j_max, time.t_1, current_time
        v_8 : float = self._segment_7_vel(time,profile_params, time.t_1)
        return v_8

    def get_velocity_at_time(self,current_time : float) -> float:
        """
        Retrieves the velocity at a specified time
        """
        return current_time
    
def main():
    profile_params : ProfileParams = ProfileParams(dist=1.0,
                                                   v_init=0.0,
                                                   v_final=10.0,
                                                   v_max=3.0,
                                                   a_max=0.8,
                                                   j_max=0.1,
                                                   delta_t=0.01)
    # type 2
    profile_params : ProfileParams = ProfileParams(10.0,
                                                   0.0,
                                                   2.0,
                                                   3,
                                                   1,
                                                   4,
                                                   0.01)


    profile : Profile = Profile(profile_params)



    total_acc : np.ndarray = np.concatenate([profile.k_1_acc,
                                        profile.k_2_acc,
                                        profile.k_3_acc,
                                        profile.k_4_acc,
                                       profile.k_5_acc,
                                       profile.k_6_acc,
                                       profile.k_7_acc,])


    total_jerk : np.ndarray = np.concatenate([profile.k_1_jerk,
                                        profile.k_2_jerk,
                                        profile.k_3_jerk,
                                        profile.k_4_jerk,
                                       profile.k_5_jerk,
                                       profile.k_6_jerk,
                                       profile.k_7_jerk,])



        # Create a figure with 3 subplots, arranged vertically
    fig, axs = plt.subplots(3, 1, figsize=(10,10))

    # Plot total_profile in the first subplot
    axs[0].plot(range(len(profile.total_profile)), profile.total_profile)
    axs[0].set_title('Total Profile')

    # Plot total_jerk in the second subplot
    axs[1].plot(range(len(total_jerk)), total_jerk)
    axs[1].set_title('Total Jerk')

    # Plot total_acc in the third subplot
    axs[2].plot(range(len(total_acc)), total_acc)
    axs[2].set_title('Total Acc')

    # Adjust the layout so that plots do not overlap
    plt.tight_layout()

    # Display the figure with the subplots
    plt.show()



if __name__ == "__main__":
    main()