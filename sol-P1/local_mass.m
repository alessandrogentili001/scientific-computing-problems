function Me = local_mass(xe)
    h = diff(xe);
    Me = h * [1/3 1/6; 1/6 1/3];
end