import s_splines

def test_profile_initial_vel():
    """
    Testing the profile if no distance is provided
    """
    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(100.0,2.0,1.0,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    assert profile.k_1[0] == 2.0, "The initial velocity is not 2.0 but: {profile.k_1[0]}"

    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(100.0,1.0,1.0,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    assert profile.k_1[0] == 1.0, "The initial velocity is not 1.0 but: {profile.k_1[0]}"

    profile_params : s_splines.ProfileParams = s_splines.ProfileParams(100.0,0.0,1.0,1,0.2,0.01)
    profile : s_splines.Profile = s_splines.Profile(profile_params)
    assert profile.k_1[0] == 0.0, "The initial velocity is not 0.0 but: {profile.k_1[0]}"