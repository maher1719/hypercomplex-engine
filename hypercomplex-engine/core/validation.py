from numbers import Integral


class Validation:
    """
    Input validation for basis elements, indices, and dimensions.
    """

    @staticmethod
    def dimension(dim) -> int:
        """
        Validates an algebra dimension exponent.

        dim must be a non-negative integer.
        dim = n means algebra dimension 2^n.
        """
        if isinstance(dim, bool) or not isinstance(dim, Integral):
            raise TypeError(f"dim must be an integer, got {type(dim).__name__}")

        dim = int(dim)

        if dim < 0:
            raise ValueError(f"dim must be >= 0, got {dim}")

        return dim

    @staticmethod
    def basis_tuple(
        data,
        allow_zero: bool = False,
        allow_eps: bool = False,
    ) -> bool:
        """
        Validates a basis-element tuple.

        Standard / Split:
            (sign, index)

        Dual:
            (sign, index, eps_flag)

        Parameters
        ----------
        data:
            The tuple to validate.

        allow_zero:
            If True, sign == 0 is allowed.
            This is needed for dual nilpotent zero products.

        allow_eps:
            If True, 3-tuples with an epsilon flag are allowed.
        """
        if not isinstance(data, tuple):
            raise TypeError(f"Expected tuple, got {type(data).__name__}")

        if len(data) not in (2, 3):
            raise ValueError(f"Basis tuple must have length 2 or 3, got {len(data)}")

        if len(data) == 3 and not allow_eps:
            raise ValueError("3-tuple epsilon form is not allowed here")
        

        sign = data[0]
        index = data[1]

        if isinstance(sign, bool) or not isinstance(sign, Integral):
            raise TypeError(f"sign must be an integer, got {type(sign).__name__}")

        if isinstance(index, bool) or not isinstance(index, Integral):
            raise TypeError(f"index must be an integer, got {type(index).__name__}")

        sign = int(sign)
        index = int(index)

        if index < 0:
            raise ValueError(f"index must be >= 0, got {index}")

        allowed_signs = (-1, 0, 1) if allow_zero else (-1, 1)

        if sign not in allowed_signs:
            raise ValueError(f"sign must be in {allowed_signs}, got {sign}")

        if len(data) == 3:
            eps = data[2]

            if isinstance(eps, bool) or not isinstance(eps, Integral):
                raise TypeError(f"eps must be an integer, got {type(eps).__name__}")

            eps = int(eps)

            if eps not in (0, 1):
                raise ValueError(f"eps must be 0 or 1, got {eps}")

        return True

    @staticmethod
    def index_in_range(index, dim: int) -> int:
        """
        Validates that index is inside [0, 2^dim - 1].
        """
        if isinstance(index, bool) or not isinstance(index, Integral):
            raise TypeError(f"index must be an integer, got {type(index).__name__}")

        index = int(index)
        dim = Validation.dimension(dim)

        size = 1 << dim

        if index < 0 or index >= size:
            raise ValueError(f"index must be in [0, {size - 1}] for dim={dim}, got {index}")

        return index