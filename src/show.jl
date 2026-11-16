# REPL displays

_sty(io, s, color; bold=false) = get(io, :color, false) ?
    string(bold ? "\e[1m" : "", "\e[38;5;$(color)m", s, "\e[0m") : string(s)
_key(io, k) = _sty(io, "  " * rpad(k, 12), 245)
_pfmt(p) = p < 0.001 ? "<0.001" : @sprintf("%.3f", p)
_stars(p) = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : p < 0.1 ? "." : ""

Base.show(io::IO, bd::BetaDisper) =
    print(io, "BetaDisper($(length(bd.distances)) observations, $(length(bd.levels)) groups)")

function Base.show(io::IO, ::MIME"text/plain", bd::BetaDisper)
    println(io, _sty(io, "BetaDisper", 117; bold=true), _sty(io, " · ", 240),
            "$(length(bd.distances)) observations in $(length(bd.levels)) groups")
    println(io, _key(io, "centre"), bd.type === :median ? "spatial median" : "centroid",
            bd.bias_adjust ? " (bias adjusted)" : "")
    nneg = count(<(0), bd.eig)
    println(io, _key(io, "axes"), "$(length(bd.eig)) principal coordinates",
            nneg > 0 ? " ($nneg with negative eigenvalues)" : "")
    println(io, _sty(io, "  mean distance to centre", 245))
    w = maximum(length ∘ string, bd.levels)
    hi = maximum(bd.group_distances)
    for (l, d) in zip(bd.levels, bd.group_distances)
        bar = repeat('━', max(1, round(Int, 20 * d / hi)))
        println(io, "    ", rpad(string(l), w), "  ", @sprintf("%.4f", d), "  ", _sty(io, bar, 117))
    end
    print(io, _sty(io, "  test with permutest(bd; pairwise=true)", 242))
end

Base.show(io::IO, t::DispersionTest) = print(io, "DispersionTest(F = ", @sprintf("%.3f", t.F), ", p = ", _pfmt(t.p), ")")

function Base.show(io::IO, ::MIME"text/plain", t::DispersionTest)
    kind = t.permutations > 0 ? "permutation test, $(t.permutations) permutations" : "F test"
    println(io, _sty(io, "DispersionTest", 117; bold=true), _sty(io, " · ", 240),
            "homogeneity of dispersions ($kind)")
    hdr = "              Df     SumSq    MeanSq        F   Pr(>F)"
    println(io, _sty(io, hdr, 245))
    st = _stars(t.p)
    println(io, "  Groups  ", lpad(t.df[1], 6), lpad(@sprintf("%.4f", t.sumsq[1]), 10),
            lpad(@sprintf("%.4f", t.meansq[1]), 10), lpad(@sprintf("%.3f", t.F), 9), lpad(_pfmt(t.p), 9), " ",
            _sty(io, st, st == "" ? 245 : 221))
    print(io, "  Residual", lpad(t.df[2], 6 - 0), lpad(@sprintf("%.4f", t.sumsq[2]), 10),
          lpad(@sprintf("%.4f", t.meansq[2]), 10))
    pw = t.pairwise
    pw === nothing && return
    println(io)
    println(io, _sty(io, "  pairwise (two-sided permutation t-tests, p adjusted: $(t.adjust))", 245))
    w = maximum(length(a) + length(b) + 3 for (a, b) in zip(pw.group1, pw.group2))
    print(io, _sty(io, "    " * rpad("pair", w) * "       t        p    p.adj", 245))
    for i in eachindex(pw.t)
        name = pw.group1[i] * " – " * pw.group2[i]
        st = _stars(pw.p_adjusted[i])
        print(io, "\n    ", rpad(name, w), lpad(@sprintf("%.3f", pw.t[i]), 8), lpad(_pfmt(pw.p[i]), 9),
              lpad(_pfmt(pw.p_adjusted[i]), 9), " ", _sty(io, st, st == "" ? 245 : 221))
    end
end
