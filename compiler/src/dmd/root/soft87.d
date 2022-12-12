module dmd.root.soft87;

/* 
 * Emulated Intel extended real floating point
 *
 *
 * Copyright: Copyright D Language Foundation 2022.
 * License:   $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
 * Authors:   Max Haughton
 */
import dmd.common.int128;
import core.stdc.stdio;


struct Real
{
    enum totalbits = 80;
    enum mantissabits = 64;
    enum expbits = 15;
    enum expmax = 2 ^^ 15 - 1;
    enum expbias = 2 ^^ (expbits - 1) - 1;
    bool sign;
    ushort exponent;
    ulong mantissa;
    void print() const
    {
        import core.stdc.stdio;

        printf("Sign: %d\n", this.sign);
        printf("Exponent: %hd (%hd)\n", this.exponent, this.exponent - expbias);
        printf("Mantissa: %zu\n", this.mantissa);
        version (all)
        {
            import std.stdio;

            writefln("%016b", this.exponent);
            writefln("%064b", this.mantissa);
            writeln(this.toNative());
        }
    }

    void printByte(ubyte x)
    {
        import core.stdc.stdio;

        for (int j = 7; j >= 0; j--)
        {
            byte bit = (x >> j) & 1;
            printf("%u", bit);
        }
    }

    @nogc:
    nothrow:
    pure:

    this(bool sign, ushort exponent, ulong mantissa) pure nothrow
    {
        this.sign = !!sign;
        this.exponent = exponent;
        this.mantissa = mantissa;
    }


    this(real x) pure nothrow
    {
        ubyte[] realBytes = (cast(ubyte*)&x)[0 .. 10];
        assert(realBytes.length == 10);
        const short front = (cast(short[]) realBytes)[$ - 1];
        import std.stdio;

        // Sign is in the high bit
        this.sign = (front >> 15) & 1;
        this.exponent = cast(ushort)(front & ~(1 << 15));
        this.mantissa = *(cast(ulong*)&realBytes[0 .. 8][0]);
    }

    real toNative() const
    {
        real x;
        import std.stdio;

        ubyte[] realBytes = (cast(ubyte*)&x)[0 .. 10];
        ushort* front = &(cast(ushort[]) realBytes)[$ - 1];
        ulong* nativeMant = (cast(ulong*)&realBytes[0 .. 8][0]);
        /+
            
        +/
        assert (!(this.exponent >> 15));
        *front = this.exponent;
        *front |= this.sign << 15;
        *nativeMant = this.mantissa;
        return x;
    }

    static typeof(this) Zero(bool sign = 0)
    {
        return typeof(this)(sign, 0, 0L);
    }

    static typeof(this) Inf(bool sign = 0)
    {
        return typeof(this)(sign, this.expmax, 0);
    }
    /// QUIET NaN
    enum QNaN = Real(0, expmax, 1UL << 63);

}
immutable(ubyte[]) x87ToBytes(const real x)
{
    ubyte* ptr = cast(ubyte*) &x;
    return ptr[0..10].idup;
}
unittest
{
    struct Should
    {
        real x;
        real y;
        bool res;
    }
    Should[] buf = [
        Should(1.0L, 1.0L, true),
        Should(1.0L, 0.9999999L, false)
    ];
    foreach(idx, sh; buf)
    {
        assert((x87ToBytes(sh.x)[] == x87ToBytes(sh.y)) == sh.res);
    }

}
void byteDiffPrint(const ubyte[] l, const ubyte[] r)
in
{
    assert(l.length == r.length);
}
do
{
    foreach(idx; 0..l.length)
    {
        import std.stdio;
        const lByte = l[idx];
        const rByte = r[idx];
        if (lByte is rByte)
        {
            //writef("%08b", );
        }
        else
        {
            debug writefln("Diff: %u", idx);
            debug writefln("%08b", lByte);
            debug writefln("%08b", rByte);
        }
        //write("\n");

    }
}

bool isNormal(const Real x)
{
    with (x)
    {
        const expRange = exponent > 0 && exponent < expmax;
        return expRange && (mantissa >> mantissabits - 1);
    }

}

unittest
{
    foreach (val; [3.14L, 1.0L])
    {
        assert(isNormal(Real(val)));
    }
}

bool isSubnormal(const Real x) pure nothrow @nogc
{
    return x.exponent == 0 && !(x.mantissa >> (x.mantissabits - 1)) && x.mantissa != 0;
}

unittest
{
    assert(!isSubnormal(Real(real.min_normal)));
    assert(isSubnormal(Real(real.min_normal / 2.0L)));
}

Real normalizeSubnormal(const Real x) pure nothrow @nogc
in
{
    assert(isSubnormal(x));
}
out (ret)
{
    debug{
        if (!isNormal(ret)) {
            x.print();
            puts("Failed!");
            ret.print();
            assert(0);
        }
    }
}
do
{
    import core.bitop : bsr;
    const shiftBy = 63 - bsr(x.mantissa);
    import std.stdio;

    return Real(x.sign, cast(short) (shiftBy), x.mantissa << shiftBy);
}

unittest
{
    puts("____________");
    import std.stdio;

    const input = Real(0, 0, 1UL << 62);
    writeln(input.toNative());
    assert(isSubnormal(input));
    input.print();
    const Real res = normalizeSubnormal(input);
    assert(res.isNormal());
    res.print();
    writeln(res.toNative());
    puts("____________");
}

Real normalizeIfSubnormal(const Real input) pure nothrow @nogc
{
    return input.isSubnormal() ? input.normalizeSubnormal : input;
}

bool isNaN(Real x) pure nothrow @nogc
{
    return x.exponent == Real.expmax && (x.mantissa >> (x.mantissabits - 2) & 0b11);
}

unittest
{
    assert(isNaN(Real(real.nan)));
    assert(isNaN(Real.QNaN));
}

bool isInf(Real x) pure nothrow
{
    return x.exponent == Real.expmax && (x.mantissa >> (x.mantissabits - 2) & 0b10);
}

unittest
{
    const real theInf = real.infinity;
    const soft = Real(theInf);
    assert(soft.toNative() is theInf);
    assert(isInf(Real(real.infinity)));
    assert(isInf(Real(-real.infinity)));
}

Real multiply(Real l, Real r) nothrow pure @nogc
{
    // debug l.print();
    // debug r.print();
    import dmd.common.int128;

    const bool outputSign = l.sign ^ r.sign;
    // Handle operations that don't require arithmetic (infs, NaN etc.)
    if (l.exponent == Real.expmax || r.exponent == Real.expmax)
    {
        // Theres a NaN
        if (isNaN(l) || isNaN(r))
        {
            // TODO: signaling nans - no one uses them but they are in the spec still and the spec is what goes.
            return Real.QNaN;
        }
        else
        {
            // Deal with infinities
            if ((l.exponent == Real.expmax && (r.exponent == 0 && r.mantissa == 0))
                || (r.exponent == Real.expmax) && (l.exponent == 0 && l.mantissa == 0))
            {
                return Real.QNaN;
            }
            else
            {
                // ( inf * -inf ) == -inf. Needs an XOR
                return Real(l.sign ^ r.sign, Real.expmax, 0L);
            }
        }
    }
    // Now for subnormals
    if (l.exponent == 0)
    {
        if (l.mantissa == 0)
        {
            return Real.Zero(outputSign);
        }
    }
    if (r.exponent == 0)
    {
        if (r.mantissa == 0)
        {
            return Real.Zero(outputSign);
        }
    }
    static auto doMul(const Real l, const Real r)
    {
        const outputSign = l.sign ^ r.sign;
        // Add then rebias exponent
        auto newExponent = l.exponent + r.exponent - Real.expbias + 1;
        Cent newMantissa = mul(Cent(l.mantissa, 0), Cent(r.mantissa, 0));
        const Cent save = newMantissa;
        import core.bitop;

        /*
            The 8087 format has a 64 bit mantissa so we don't have to do anything clever
            in or before the integer multiply.
        */
       /*  while (!(newMantissa.hi >> 63))
        {
            import std.stdio;

            writefln("-%064b %064b", newMantissa.hi, newMantissa.lo);
            newMantissa = shl1(newMantissa);
            writefln("+%064b %064b", newMantissa.hi, newMantissa.lo);
            --newExponent;
        } */

        import std.stdio;

        const shiftBy = 63 - bsr(save.hi);
        const tryShift = shl(save, shiftBy);
        // writefln("%064b %064b", save.hi, save.lo);
        // writefln("%064b %064b", tryShift.hi, tryShift.lo);
        /*
            The aim here is to support round to nearest
        */
        const guardBit = (tryShift.lo & (1UL << (Real.mantissabits - 1))) != 0;
        const stickyBit = (tryShift.lo & ((1UL << (Real.mantissabits - 1)) - 1)) != 0;

        const ulong finalMantissa = tryShift.hi + (guardBit && (stickyBit || (tryShift.hi & 1)));
        const finalExp = newExponent > shiftBy ? newExponent - shiftBy : 0;
        // writefln("Exp is %u", finalExp);

        return Real(outputSign, cast(ushort) finalExp, finalMantissa);
    }

    return doMul(normalizeIfSubnormal(l), normalizeIfSubnormal(r));
}

unittest
{
    assert(multiply(Real.Zero(), Real.Inf()) is Real.QNaN);
    assert(multiply(Real.Inf(), Real.Zero()) is Real.QNaN);
    assert(multiply(Real.Zero(), Real.Zero()) is Real.Zero());
}

unittest
{
    real[2][] list;
    list = [
        [1.0L, 1.0L],
        [400.0L, 0.0001L],
        [1 / 3.0L, 10000.0L],
        [1.0 / 3.0L, 1.0 / 11.0L],
        [real.min_normal / 2, real.min_normal],
        [real.max, real.min_normal],
        [314.00L, 1.000L]
    ];
    foreach (arr; list)
    {
        const real x = arr[0];
        const real y = arr[1];
        auto rb = Real(x);
        const o = rb.toNative();
        auto bb = Real(y);
        {
            const xx = multiply(rb, bb);
            const real res = xx.toNative();
            const real native = x * y;
            import std.stdio;
            writefln("%f, %f", res, native);
            const resBytes = x87ToBytes(res);
            const nativeBytes = x87ToBytes(res);
            if (!(resBytes[] == nativeBytes[]))
            {
                import std.stdio;
                writeln([x, y]);
                xx.print();

                writeln(res);
                writeln(native);
                byteDiffPrint(x87ToBytes(res), x87ToBytes(native));
            }
        }
    }
}
