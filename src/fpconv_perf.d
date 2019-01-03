import fpconv_ctfe;
import vmm;

import core.stdc.stdio;

extern (C) void usleep(uint usecs);

double ticks_to_usecs = 0.0;

char[] itos(uint n)
{
    char[] result = [];
    const len = 10;
    result.length = len;
    uint i = len - 1;

    while (n > 10)
    {
        result[i--] = cast(char)('0' + (n % 10));
        n /= 10;
    }
    result[i] = cast(char)('0' + (n % 10));

    return result[i .. $];
}

char[] TRACEPOINT_ID_STRING(char[] name = null,
        char[] file = cast(char[]) __FILE__, uint line = __LINE__)
{
    if (name)
    {
        return file ~ "_L" ~ itos(line) ~ "_" ~ name;
    }
    else
    {
        return file ~ "_L" ~ itos(line) ~ "_";
    }
}

/*
enum proile_string = 
q{enum __bgn_counter__ = mixin("uint 
 = void; asm { rdtsc; mov EAX");
};
*/

void init_ticks_to_wallclock()
{
LBegin:
    static long begin = void;
    begin = rdtsc();

    usleep(59817);

    static long end = void;
    end = rdtsc();

    double real_elapsed_approx = (59817 - 20);

    if (end < begin)
    {
        // we wrapped around. repeat :)
        goto LBegin;
    }

    ticks_to_usecs = real_elapsed_approx / (cast(double)(end - begin));

    printf("rdtsc currently at about : %f mhz \n".ptr, 1 / ticks_to_usecs);
    printf("precision of uint will last for %f seconds \n".ptr,
            ((uint.max / 1_000_000_000) * ticks_to_usecs));

}
/+
extern (C) static uint rdtsc()
{
    asm { naked; rdtsc; ret; }
}
+/

ubyte[] GetPageBytes(int n)
{
    ubyte[] result;
    result = cast(ubyte[]) allocate(null, n);
    commit(result);
    return result;
}

void SetPageExecutable(void[] page)
{
    protect(page, MemoryAccess.execute);
}



extern (C) long function() rdtsc;
void main()
{

    // assert(double(.670065995362982058) == double(670065995362982058));
    // Get one Page of ExecutableMemeory
    ubyte[] execmem = GetPageBytes(4096);

    double[] doubles =
    [
        0.6700659953629821, 0.8618356877334371, 0.8424251474360873,
        0.32723744881029593, 0.1753466728325745, 0.3257012747075798,
        0.32081499588437645, 0.6719597837196847, 0.48548373558652463,
        0.5360050465082955, 0.9106717974818522, 0.2673529950328973, 
        0.7259865338863339, 0.4262670635969924, 0.9456009305779313,
        0.3111753234314665
    ];

    execmem[0 .. 3] = [0x0f, 0x31, 0xc3];
    SetPageExecutable(execmem);

    rdtsc = cast(typeof(rdtsc))&execmem[0];

    init_ticks_to_wallclock();

    char[24][2048] results;
    // char[24][128] dota_results;
    int[results.length] resultLengths;


        foreach (i; 0 .. results.length)
        {
            double d = doubles[i % doubles.length];
        }



    foreach (i; 0 .. results.length)
    {
        printf("%.*s\n", resultLengths[i], &results[i][0]);
    }

    //{
        auto dtoa_begin_ticks = GetTicks();

        foreach (i; 0 .. results.length)
        {
            double d = doubles[i % doubles.length];
            resultLengths[i] = fpconv_dtoa(d, results[i]);
        }

        auto dtoa_end_ticks = GetTicks();
    //}



    foreach (i; 0 .. results.length)
    {
        // printf("%.*s\n", resultLengths[i], &results[i][0]);
    }


    // {
        auto snpf_begin_ticks = GetTicks();

        foreach (i; 0 .. results.length)
        {
            double d = doubles[i % doubles.length];
            resultLengths[i] = snprintf(&results[i][0], 24, "%.24f", d);
        }

        auto snpf_end_ticks = GetTicks();
   // }

    foreach (i; 0 .. results.length)
    {
        // printf("%.*s\n", resultLengths[i], &results[i][0]);
    }


    //{
        auto stb_snpf_begin_ticks = GetTicks();

        foreach (i; 0 .. results.length)
        {
            double d = doubles[i % doubles.length];
            const string s = double2str(d, 18);
            const l = (s.length < 24 ? cast(byte) s.length : byte(24));
            results[i][0 .. l] = s[0 .. l];
            resultLengths[i] = l;
        }

        auto stb_snpf_end_ticks = GetTicks();
    //}


    printf("\n\nConverting %d floats with dota took %f us\n",
            results.length * 1, (cast(double)dtoa_end_ticks - dtoa_begin_ticks) * ticks_to_usecs);

    printf("\n\nConverting %d floats with snpf took %f us\n",
            results.length, (cast(double)snpf_end_ticks - snpf_begin_ticks) * ticks_to_usecs);

    printf("\n\nConverting %d floats with stb_snpf took %f us\n",
            results.length, (cast(double)stb_snpf_end_ticks - stb_snpf_begin_ticks) * ticks_to_usecs);

}

alias GetTicks = rdtsc;

/// generates a uniform random double
double splitmix64d(double x) pure
{
    ulong z = *cast(ulong*)&x;
    z += 0x9e3779b97f4a7c15;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    z = z ^ (z >> 31);
    return z * 5.42101e-20;
}

struct urd
{
    int n;
    double state;

    double front()
    {
        return state;
    }

    void popFront()
    {
        state = splitmix64d(state);
    }

    bool empty()
    {
        return n-- == 0;
    }
}


///////////////////////////////////////////// ketmars stuff ////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////// //
public string double2str (double d, int fracdigs=-1) nothrow @trusted {
  char[STBSP__NUMSZ] num;
  uint fdc = (fracdigs >= 0 ? cast(uint)fracdigs : 6);
  uint olen;
  int dpos;
  const(char)* start;
  int sign = stbsp__real_to_str(&start, &olen, num.ptr, &dpos, d, fdc);
  if (dpos == STBSP__SPECIAL) {
    return start[0..olen].idup;
  } else {
    while (olen > dpos && start[olen-1] == '0') --olen;
    uint reslen = (sign ? 1 : 0)+olen+(dpos == 0 ? 2 : 1);
    char[] res = new char[](reslen+(olen-dpos < fdc ? fdc-(olen-dpos) : 0));
    res[] = '0'; // why not?
    uint rpos = 0;
    if (sign) res.ptr[rpos++] = '-';
    if (dpos > 0) { res.ptr[rpos..rpos+dpos] = start[0..dpos]; rpos += dpos; }
    if (fdc != 0) {
      res.ptr[rpos++] = '.';
      res.ptr[rpos..rpos+(olen-dpos)] = start[dpos..olen];
      rpos += olen-dpos;
      if (fracdigs > 0 && olen-dpos < fdc) rpos += fdc-(olen-dpos);
    }
    assert(rpos <= res.length);
    return cast(string)res.ptr[0..rpos]; // it is safe to cast here
  }
}


// ////////////////////////////////////////////////////////////////////////// //
private:

public enum STBSP__SPECIAL = 0x7000;

public enum STBSP__NUMSZ = 512; // big enough for e308 (with commas) or e-307
//char[STBSP__NUMSZ] num;


// ////////////////////////////////////////////////////////////////////////// //
static immutable double[23] stbsp__bot = [
  0x1p+0,0x1.4p+3,0x1.9p+6,0x1.f4p+9,0x1.388p+13,0x1.86ap+16,0x1.e848p+19,0x1.312dp+23,
  0x1.7d784p+26,0x1.dcd65p+29,0x1.2a05f2p+33,0x1.74876e8p+36,0x1.d1a94a2p+39,0x1.2309ce54p+43,
  0x1.6bcc41e9p+46,0x1.c6bf52634p+49,0x1.1c37937e08p+53,0x1.6345785d8ap+56,0x1.bc16d674ec8p+59,
  0x1.158e460913dp+63,0x1.5af1d78b58c4p+66,0x1.b1ae4d6e2ef5p+69,0x1.0f0cf064dd592p+73,
];

static immutable double[22] stbsp__negbot = [
  0x1.999999999999ap-4,0x1.47ae147ae147bp-7,0x1.0624dd2f1a9fcp-10,0x1.a36e2eb1c432dp-14,
  0x1.4f8b588e368f1p-17,0x1.0c6f7a0b5ed8dp-20,0x1.ad7f29abcaf48p-24,0x1.5798ee2308c3ap-27,
  0x1.12e0be826d695p-30,0x1.b7cdfd9d7bdbbp-34,0x1.5fd7fe1796495p-37,0x1.19799812dea11p-40,
  0x1.c25c268497682p-44,0x1.6849b86a12b9bp-47,0x1.203af9ee75616p-50,0x1.cd2b297d889bcp-54,
  0x1.70ef54646d497p-57,0x1.2725dd1d243acp-60,0x1.d83c94fb6d2acp-64,0x1.79ca10c924223p-67,
  0x1.2e3b40a0e9b4fp-70,0x1.e392010175ee6p-74,
];

static immutable double[22] stbsp__negboterr = [
  -0x1.999999999999ap-58,-0x1.eb851eb851eb8p-63,-0x1.89374bc6a7efap-66,-0x1.6a161e4f765fep-68,
  -0x1.ee78183f91e64p-71,0x1.b5a63f9a49c2cp-75,0x1.5e1e99483b023p-78,-0x1.03023df2d4c94p-82,
  -0x1.34674bfabb83bp-84,-0x1.20a5465df8d2cp-88,0x1.7f7bc7b4d28aap-91,0x1.97f27f0f6e886p-96,
  -0x1.ecd79a5a0df95p-99,0x1.ea70909833de7p-107,-0x1.937831647f5ap-104,0x1.5b4c2ebe68799p-109,
  -0x1.db7b2080a3029p-111,-0x1.7c628066e8ceep-114,0x1.a52b31e9e3d07p-119,0x1.75447a5d8e536p-121,
  0x1.f769fb7e0b75ep-124,-0x1.a7566d9cba769p-128,
];

static immutable double[13] stbsp__top = [
  0x1.52d02c7e14af6p+76,0x1.c06a5ec5433c6p+152,0x1.28bc8abe49f64p+229,0x1.88ba3bf284e24p+305,
  0x1.03e29f5c2b18cp+382,0x1.57f48bb41db7cp+458,0x1.c73892ecbfbf4p+534,0x1.2d3d6f88f0b3dp+611,
  0x1.8eb0138858d0ap+687,0x1.07d457124123dp+764,0x1.5d2ce55747a18p+840,0x1.ce2137f743382p+916,
  0x1.31cfd3999f7bp+993,
];

static immutable double[13] stbsp__negtop = [
  0x1.82db34012b251p-77,0x1.244ce242c5561p-153,0x1.b9b6364f30304p-230,0x1.4dbf7b3f71cb7p-306,
  0x1.f8587e7083e3p-383,0x1.7d12a4670c123p-459,0x1.1fee341fc585dp-535,0x1.b31bb5dc320d2p-612,
  0x1.48c22ca71a1bdp-688,0x1.f0ce4839198dbp-765,0x1.77603725064a8p-841,0x1.1ba03f5b21p-917,
  0x1.ac9a7b3b7302fp-994,
];

static immutable double[13] stbsp__toperr = [
  0x1p+23,0x1.bb542c80deb4p+95,-0x1.83b80b9aab60ap+175,-0x1.32e22d17a166cp+251,
  -0x1.23606902e180ep+326,-0x1.96fb782462e87p+403,-0x1.358952c0bd011p+480,-0x1.78c1376a34b6cp+555,
  -0x1.17569fc243adfp+633,-0x1.d9365a897aaa6p+710,0x1.9050c256123ap+786,-0x1.b1799d76cc7a6p+862,
  -0x1.213fe39571a38p+939,
];

static immutable double[13] stbsp__negtoperr = [
  0x1.13badb829e079p-131,-0x1.e46a98d3d9f64p-209,0x1.227c7218a2b65p-284,0x1.1d96999aa01e9p-362,
  -0x1.cc2229efc3962p-437,-0x1.cd04a2263407ap-513,-0x1.23b80f187a157p-590,-0x1.c4e22914ed912p-666,
  0x1.bc296cdf42f82p-742,-0x1.9f9e7f4e16fe1p-819,-0x1.aeb0a72a8902ap-895,-0x1.e228e12c13408p-971,
  0x0.0000000fa1259p-1022,
];


static immutable ulong[20] stbsp__powten = [
   1UL,
   10UL,
   100UL,
   1000UL,
   10000UL,
   100000UL,
   1000000UL,
   10000000UL,
   100000000UL,
   1000000000UL,
   10000000000UL,
   100000000000UL,
   1000000000000UL,
   10000000000000UL,
   100000000000000UL,
   1000000000000000UL,
   10000000000000000UL,
   100000000000000000UL,
   1000000000000000000UL,
   10000000000000000000UL
];
enum stbsp__tento19th = 1000000000000000000UL;

static immutable string stbsp__digitpair = "00010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899";


// get float info; returns `-1` for negatives (including negative zero), and `0` for positives
public int stbsp__real_to_parts (long* bits, int* expo, double value) nothrow @trusted @nogc {
  // load value and round at the frac_digits
  double d = value;
  long b = *(cast(const(long)*)&d);
  //STBSP__COPYFP(b, d);

  *bits = b&(((cast(ulong)1)<<52)-1);
  *expo = cast(int)(((b>>52)&2047)-1023);

  return cast(int)(b>>63);
}


void stbsp__ddmulthi (ref double oh, ref double ol, in double xh, in double yh) nothrow @trusted @nogc {
  double ahi = 0, alo, bhi = 0, blo;
  long bt;
  oh = xh * yh;
  //STBSP__COPYFP(bt, xh);
  bt = *(cast(const(long)*)&xh);
  bt &= ((~cast(ulong)0) << 27);
  //STBSP__COPYFP(ahi, bt);
  ahi = *(cast(const(double)*)&bt);
  alo = xh - ahi;
  //STBSP__COPYFP(bt, yh);
  bt = *(cast(const(long)*)&yh);
  bt &= ((~cast(ulong)0) << 27);
  //STBSP__COPYFP(bhi, bt);
  bhi = *(cast(const(double)*)&bt);
  blo = yh - bhi;
  ol = ((ahi * bhi - oh) + ahi * blo + alo * bhi) + alo * blo;
}

void stbsp__ddtoS64 (ref long ob, in double xh, in double xl, in double ph) nothrow @trusted @nogc {
  double ahi = 0, alo, vh, t;
  ob = cast(long)ph;
  vh = cast(double)ob;
  ahi = (xh - vh);
  t = (ahi - xh);
  alo = (xh - (ahi - t)) - (vh + t);
  ob += cast(long)(ahi + alo + xl);
}

void stbsp__ddrenorm (ref double oh, ref double ol) nothrow @trusted @nogc {
  double s = oh + ol;
  ol = ol - (s - oh);
  oh = s;
}

void stbsp__ddmultlo (in double oh, ref double ol, in double xh, in double xl, in double yh, in double yl) nothrow @trusted @nogc { ol = ol + (xh * yl + xl * yh); }

void stbsp__ddmultlos (in double oh, ref double ol, in double xh, in double yl) nothrow @trusted @nogc { ol = ol + (xh * yl); }

// power can be -323 to +350
void stbsp__raise_to_power10 (double* ohi, double* olo, double d, int power) nothrow @trusted @nogc {
  double ph, pl;
  if (power >= 0 && power <= 22) {
    stbsp__ddmulthi(ph, pl, d, stbsp__bot[power]);
  } else {
    int e, et, eb;
    double p2h, p2l;

    e = power;
    if (power < 0) e = -e;
    et = (e * 0x2c9) >> 14; /* %23 */
    if (et > 13) et = 13;
    eb = e - (et * 23);

    ph = d;
    pl = 0.0;
    if (power < 0) {
      if (eb) {
        --eb;
        stbsp__ddmulthi(ph, pl, d, stbsp__negbot[eb]);
        stbsp__ddmultlos(ph, pl, d, stbsp__negboterr[eb]);
      }
      if (et) {
        stbsp__ddrenorm(ph, pl);
        --et;
        stbsp__ddmulthi(p2h, p2l, ph, stbsp__negtop[et]);
        stbsp__ddmultlo(p2h, p2l, ph, pl, stbsp__negtop[et], stbsp__negtoperr[et]);
        ph = p2h;
        pl = p2l;
      }
    } else {
      if (eb) {
        e = eb;
        if (eb > 22) eb = 22;
        e -= eb;
        stbsp__ddmulthi(ph, pl, d, stbsp__bot[eb]);
        if (e) {
          stbsp__ddrenorm(ph, pl);
          stbsp__ddmulthi(p2h, p2l, ph, stbsp__bot[e]);
          stbsp__ddmultlos(p2h, p2l, stbsp__bot[e], pl);
          ph = p2h;
          pl = p2l;
        }
      }
      if (et) {
        stbsp__ddrenorm(ph, pl);
        --et;
        stbsp__ddmulthi(p2h, p2l, ph, stbsp__top[et]);
        stbsp__ddmultlo(p2h, p2l, ph, pl, stbsp__top[et], stbsp__toperr[et]);
        ph = p2h;
        pl = p2l;
      }
    }
  }
  stbsp__ddrenorm(ph, pl);
  *ohi = ph;
  *olo = pl;
}


// given a float value, returns the significant bits in bits, and the position of the decimal point in decimal_pos.
// +INF/-INF and NAN are specified by special values returned in the decimal_pos parameter.
// frac_digits is absolute normally, but if you want from first significant digits (got %g and %e), or in 0x80000000.
// returns `-1` for negatives (including negative zero), and `0` for positives.
public int stbsp__real_to_str (const(char)** start, uint* len, char* outbuf, int* decimal_pos, double value, uint frac_digits=6) nothrow @trusted @nogc {
  double d;
  long bits = 0;
  int expo, e, ng, tens;

  d = value;
  //STBSP__COPYFP(bits, d);
  bits = *(cast(const(long)*)&d);
  expo = cast(int)((bits >> 52) & 2047);
  ng = cast(int)(bits >> 63);
  if (ng) d = -d;

  // is nan or inf?
  if (expo == 2047) {
    *start = (bits & (((cast(ulong)1) << 52) - 1)) ? "NaN" : "Inf";
    *decimal_pos = STBSP__SPECIAL;
    *len = 3;
    return ng;
  }

  // is zero or denormal?
  if (expo == 0) {
    if ((bits << 1) == 0) {
      // do zero
      *decimal_pos = 1;
      *start = outbuf;
      outbuf[0] = '0';
      *len = 1;
      return ng;
    }
    // find the right expo for denormals
    {
      long v = (cast(ulong)1) << 51;
      while ((bits & v) == 0) {
        --expo;
        v >>= 1;
      }
    }
  }

  // find the decimal exponent as well as the decimal bits of the value
  {
    double ph, pl;

    // log10 estimate - very specifically tweaked to hit or undershoot by no more than 1 of log10 of all expos 1..2046
    tens = expo - 1023;
    tens = (tens < 0) ? ((tens * 617) / 2048) : (((tens * 1233) / 4096) + 1);

    // move the significant bits into position and stick them into an int
    stbsp__raise_to_power10(&ph, &pl, d, 18 - tens);

    // get full as much precision from double-double as possible
    stbsp__ddtoS64(bits, ph, pl, ph);

    // check if we undershot
    if ((cast(ulong)bits) >= stbsp__tento19th) ++tens;
  }

  // now do the rounding in integer land
  frac_digits = (frac_digits & 0x80000000) ? ((frac_digits & 0x7ffffff) + 1) : (tens + frac_digits);
  if ((frac_digits < 24)) {
    uint dg = 1;
    if (cast(ulong)bits >= stbsp__powten[9]) dg = 10;
    while (cast(ulong)bits >= stbsp__powten[dg]) {
      ++dg;
      if (dg == 20) goto noround;
    }
    if (frac_digits < dg) {
      ulong r;
      // add 0.5 at the right position and round
      e = dg - frac_digits;
      if (cast(uint)e >= 24) goto noround;
      r = stbsp__powten[e];
      bits = bits + (r / 2);
      if (cast(ulong)bits >= stbsp__powten[dg])
         ++tens;
      bits /= r;
    }
  noround:;
  }

  // kill long trailing runs of zeros
  if (bits) {
    uint n;
    for (;;) {
      if (bits <= 0xffffffff) break;
      if (bits % 1000) goto donez;
      bits /= 1000;
    }
    n = cast(uint)bits;
    while ((n % 1000) == 0) n /= 1000;
    bits = n;
  donez:;
  }

  // convert to string
  outbuf += 64;
  e = 0;
  for (;;) {
    uint n;
    char* o = outbuf-8;
    // do the conversion in chunks of U32s (avoid most 64-bit divides, worth it, constant denomiators be damned)
    if (bits >= 100000000) {
      n = cast(uint)(bits % 100000000);
      bits /= 100000000;
    } else {
      n = cast(uint)bits;
      bits = 0;
    }
    while (n) {
      outbuf -= 2;
      if (__ctfe) {
        outbuf[0..2] = stbsp__digitpair[(n%100)*2..(n%100)*2+2];
      } else {
        *cast(ushort*)outbuf = *cast(ushort*)&stbsp__digitpair[(n%100)*2];
      }
      n /= 100;
      e += 2;
    }
    if (bits == 0) {
      if (e && outbuf[0] == '0') {
        ++outbuf;
        --e;
      }
      break;
    }
    while (outbuf != o) {
      *--outbuf = '0';
      ++e;
    }
  }

  *decimal_pos = tens;
  *start = outbuf;
  *len = e;
  return ng;
}
