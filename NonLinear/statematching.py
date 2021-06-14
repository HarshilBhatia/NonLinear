def W_forward(u, v, w, e, eta):

    """
    Computes W's next iteration.
    Parameters: 
    u (Double)
    v (Double)
    w (Double)
    e (Double)
    eta (Double)

    Returns : iterated value of w
    """

    t1_n = eta * (w - 1) ** 2
    t1_d = (
        (e ** 2 - 1) * u ** 2
        + (e ** 2 - 1) * v ** 2
        + 2 * e ** 2 * w ** 2
        + 2 * (e ** 2 - 1) * w
        + 2
    )

    t1 = t1_n / t1_d

    t2_n = eta * (
        ((e ** 2) - 1) * u ** 2
        + (e ** 2 - 1) * v ** 2
        + (w + 1) * ((2 * (e ** 2) - 1) * w + 1)
    )

    t2 = t2_n / t1_d

    t3_n = (eta - 1) * (u ** 2 + v ** 2 - w ** 2 + 1)
    t3_d = (
        (e ** 2 - 1) * u ** 2 + (e ** 2 - 1) * v ** 2 + 2 * (w + 1) * ((e ** 2) * w - 1)
    )
    t3 = t3_n / t3_d

    t4_n = (eta - 1) * (e ** 2 * (u ** 2 + v ** 2 + 2 * w * (w + 1)) - (w + 1) ** 2)
    t4 = t4_n / t3_d

    return -t1 + t2


# def W_temp(w,e):
#     t1_n = (1-e)


def U_forward(u, v, w, e, eta):
    """
    Computes u's next iteration.
    Parameters: 
    u (Double)
    v (Double)
    w (Double)
    e (Double)
    eta (Double)

    Returns : iterated value of u
    """
    t1_n = 2 * e * eta * (u ** 2 - v ** 2)
    t1_d = (
        (e ** 2 - 1) * u ** 2
        + (e ** 2 - 1) * v ** 2
        + 2 * e ** 2 * w ** 2
        + 2 * (e ** 2 - 1) * 2
        + 2
    )
    t1 = t1_n / t1_d

    t2_n = 2 * (2 - 2 * e ** 2) ** 0.5 * (eta - 1) * u * (1 + w)
    t2_d = (
        (e ** 2 - 1) * u ** 2 + (e ** 2 - 1) * v ** 2 + 2 * (w + 1) * ((e ** 2) * w - 1)
    )
    t2 = t2_n / t2_d

    return t1 + t2


def V_forward(u, v, w, e, eta):
    """
    Computes v's next iteration.
    Parameters: 
    u (Double)
    v (Double)
    w (Double)
    e (Double)
    eta (Double)

    Returns : iterated value of v
    """
    t1_n = 4 * e * eta * u * v
    t1_d = (
        (e ** 2 - 1) * u ** 2
        + (e ** 2 - 1) * v ** 2
        + 2 * e ** 2 * w ** 2
        + 2 * (e ** 2 - 1) * 2
        + 2
    )
    t1 = t1_n / t1_d

    t2_n = 2 * (2 - 2 * e ** 2) ** 0.5 * (eta - 1) * v * (1 + w)
    t2_d = (
        (e ** 2 - 1) * u ** 2 + (e ** 2 - 1) * v ** 2 + 2 * (w + 1) * ((e ** 2) * w - 1)
    )
    t2 = t2_n / t2_d

    return t1 + t2

