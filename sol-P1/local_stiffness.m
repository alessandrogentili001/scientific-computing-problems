function Ke = local_stiffness(xe)
    h = diff(xe);
    Ke = [1 -1; -1 1] / h;
end