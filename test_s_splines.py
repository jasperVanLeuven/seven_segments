import s_splines

# VELOCITY START
def test_profile_initial_vel():
    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(20.0,2.0,4.0,10,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    assert profile.k_1_vel[0] == 2.0, f"The initial velocity is not 2.0 but: {profile.k_1_vel[0]}"

def test_profile_initial_vel_1():
    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(2.0,1.0,5.0,10,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    assert profile.k_1_vel[0] == 1.0, f"The initial velocity is not 1.0 but: {profile.k_1_vel[0]}"

def test_profile_initial_vel_2():
    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(40.0,0.0,10.0,10,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    assert profile.k_1_vel[0] == 0.0, f"The initial velocity is not 0.0 but: {profile.k_1_vel[0]}"

# VELOCITY END
def test_profile_end_vel():
    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(10.0,2.0,4.0,10,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    real_final : float = profile.k_7_vel[-1]
    calc_final : float = profile_params.v_final
    assert abs(real_final - calc_final) < 0.0001, f"The final velocity is not {calc_final} but: {real_final}"

def test_profile_end_vel_1():
    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(20.0,1.0,5.0,10,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    real_final : float = profile.k_7_vel[-1]
    calc_final : float = profile_params.v_final
    assert abs(real_final - calc_final) < 0.0001, f"The final velocity is not {calc_final} but: {real_final}"

def test_profile_end_vel_2():
    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(20.0,0.0,10.0,10,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    real_final : float = profile.k_7_vel[-1]
    calc_final : float = profile_params.v_final
    assert abs(real_final - calc_final) < 0.0001, f"The final velocity is not {calc_final} but: {real_final}"