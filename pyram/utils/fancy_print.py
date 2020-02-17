from sty import fg, bg, ef, rs


def printw(s, color="orange"):
    """
    print warnings in a different color

    You can register a custom color by doing
    >>> fg.orange = Rule(Render.rgb_fg, 255, 150, 50)


    Refer https://github.com/feluxe/sty for more details.
    """
    print(fg.li_green + "[Warning]" + s + fg.rs)
