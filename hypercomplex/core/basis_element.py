class BasisElement:
    """
    Basis elements are represented as plain tuples.

    Standard / Split:
        (sign, index)

    Dual:
        (sign, index, eps_flag)

    sign:
        -1, 0, or +1
        sign == 0 means the zero element.

    index:
        basis index k, meaning e_k.

    eps_flag:
        0 or 1, used only for dual numbers.
    """

    @staticmethod
    def make(sign: int, index: int, eps: int | None = None) -> tuple:
        if eps is None:
            return (int(sign), int(index))
        return (int(sign), int(index), int(eps))

    @staticmethod
    def zero(eps: int = 0) -> tuple:
        return (0, 0, int(eps))

    @staticmethod
    def sign(element: tuple) -> int:
        return int(element[0])

    @staticmethod
    def index(element: tuple) -> int:
        return int(element[1])

    @staticmethod
    def eps(element: tuple) -> int:
        if len(element) >= 3:
            return int(element[2])
        return 0

    @staticmethod
    def is_zero(element: tuple) -> bool:
        return int(element[0]) == 0

    @staticmethod
    def negate(element: tuple) -> tuple:
        if len(element) == 2:
            return (-int(element[0]), int(element[1]))
        return (-int(element[0]), int(element[1]), int(element[2]))

    @staticmethod
    def with_eps(element: tuple, eps: int) -> tuple:
        return (int(element[0]), int(element[1]), int(eps))