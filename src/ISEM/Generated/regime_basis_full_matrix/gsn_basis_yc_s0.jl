# Compact Horner evaluation; coefficients are constructed as in the Mathematica export.
# Identical constructor expressions reuse an earlier coefficient without reassociation.
# Generated from the Mathematica export script direct_gsn_basis_full_matrix_export.wl (not distributed with the package).
# Each monomial's polynomial in c = cos(nu) carries its exact factor (1 - c)^k as u^k.
# Split scope: basis=yc, spin=0.
# Do not edit by hand.

module GSNBasis_yc_s0
using ...CoefficientTables: _evaluate_polynomial_blocks!

Base.Experimental.@compiler_options optimize=0

export construct_gsn_basis, eval_gsn_basis!, eval_gsn_basis_variable!, ordinary_value_count

const LEN_GSN_BASIS_yc_s0 = 136
const N_GSN_BASIS_yc_s0 = 43

# Shared-preamble constructor blocks: split_direct_gsn_generated_constructors.jl

@noinline function _construct_gsn_basis_yc_s0_block_1!(coeffs, cp2, cp3, cp4, u)
    coeffs[1] = ComplexF64(32*cp4)
    coeffs[2] = ComplexF64(16*cp3 - 64*cp4)
    coeffs[3] = ComplexF64(-32*cp3*u)
    coeffs[4] = ComplexF64(-8*cp2*u)
    coeffs[5] = ComplexF64(-16*cp3 + 64*cp4)
    coeffs[6] = ComplexF64(64*cp3*u)
    coeffs[7] = ComplexF64(24*cp2*u)
    coeffs[8] = coeffs[3]
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_2!(coeffs, cp2, cp3, cp4, u)
    coeffs[9] = ComplexF64(-24*cp2*u)
    coeffs[10] = ComplexF64(8*cp2*u)
    coeffs[11] = ComplexF64(0)
    coeffs[12] = ComplexF64(-16*cp4)
    coeffs[13] = ComplexF64(-16*cp3 + 32*cp4)
    coeffs[14] = ComplexF64(-8*cp2*u + 16*cp3*u)
    coeffs[15] = coeffs[10]
    coeffs[16] = ComplexF64(16*cp4)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_3!(coeffs, cp2, cp3, cp4, u)
    coeffs[17] = ComplexF64(32*cp3 - 64*cp4)
    coeffs[18] = ComplexF64(24*cp2*u - 48*cp3*u)
    coeffs[19] = ComplexF64(-32*cp2*u)
    coeffs[20] = coeffs[13]
    coeffs[21] = ComplexF64(-24*cp2*u + 48*cp3*u)
    coeffs[22] = ComplexF64(48*cp2*u)
    coeffs[23] = ComplexF64(8*cp2*u - 16*cp3*u)
    coeffs[24] = coeffs[19]
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_4!(coeffs, lambda, m, omega, m2, o2, sn, cp4, cp5, cp6, cp7, cp8, cp9, cp10, u, up2, up3)
    coeffs[25] = coeffs[10]
    coeffs[26] = ComplexF64(1024*o2*cp10)
    coeffs[27] = ComplexF64(4096*o2*cp9*u)
    coeffs[28] = ComplexF64(8192*o2*cp8*u - 6144*o2*cp9*u - 256*lambda*cp8 - 512*sn*m*omega*cp8)
    coeffs[29] = ComplexF64(-256*cp7 + 10240*o2*cp7*up2 - 4096*o2*cp8*up2 - 512*lambda*cp7 + 768*lambda*cp8 - 1536*sn*m*omega*cp7*u)
    coeffs[30] = ComplexF64(-448*cp6 + 640*cp7 + 64*cp8 + 8704*o2*cp6*up2 - 8192*o2*cp7*up2 + 1024*o2*cp8*up2 + 64*m2*cp6*u + 64*m2*cp7*u - 512*lambda*cp6*u + 768*lambda*cp7*u - 2304*sn*m*omega*cp6*u + 1536*sn*m*omega*cp7*u)
    coeffs[31] = ComplexF64(-128*cp5*u + 768*cp6*u + 128*cp7*u + 5120*o2*cp5*up3 - 2048*o2*cp6*up3 + 128*m2*cp5*up2 + 128*m2*cp6*up2 - 256*lambda*cp5*u + 768*lambda*cp6*u - 256*lambda*cp7*u - 2048*sn*m*omega*cp5*up2 + 512*sn*m*omega*cp6*up2)
    coeffs[32] = ComplexF64(192*cp5*u - 384*cp6*u - 64*cp7*u + 2048*o2*cp4*up3 - 1536*o2*cp5*up3 + 128*m2*cp4*up2 + 64*m2*cp5*up2 - 64*m2*cp6*up2 - 64*lambda*cp4*up2 + 256*lambda*cp5*up2 - 1152*sn*m*omega*cp4*up2 + 768*sn*m*omega*cp5*up2)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_5!(coeffs, lambda, m, omega, m2, o2, sn, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9, u, up2, up3, up4)
    coeffs[33] = ComplexF64(-64*cp5*up2 + 512*o2*cp3*up4 + 64*m2*cp3*up3 + 64*m2*cp4*up3 + 64*lambda*cp4*up2 - 384*sn*m*omega*cp3*up3)
    coeffs[34] = ComplexF64(64*o2*cp2*up4 + 16*m2*cp2*up3 + 16*m2*cp3*up3 - 64*sn*m*omega*cp2*up3)
    coeffs[35] = ComplexF64(-4096*o2*cp9*u)
    coeffs[36] = ComplexF64(-16384*o2*cp8*u + 12288*o2*cp9*u + 512*lambda*cp8 + 1024*sn*m*omega*cp8)
    coeffs[37] = ComplexF64(768*cp7 - 30720*o2*cp7*up2 + 12288*o2*cp8*up2 + 1536*lambda*cp7 - 2304*lambda*cp8 + 4608*sn*m*omega*cp7*u)
    coeffs[38] = ComplexF64(1792*cp6 - 2560*cp7 - 256*cp8 - 34816*o2*cp6*up2 + 32768*o2*cp7*up2 - 4096*o2*cp8*up2 - 256*m2*cp6*u - 256*m2*cp7*u + 2048*lambda*cp6*u - 3072*lambda*cp7*u + 9216*sn*m*omega*cp6*u - 6144*sn*m*omega*cp7*u)
    coeffs[39] = ComplexF64(640*cp5*u - 3840*cp6*u - 640*cp7*u - 25600*o2*cp5*up3 + 10240*o2*cp6*up3 - 640*m2*cp5*up2 - 640*m2*cp6*up2 + 1280*lambda*cp5*u - 3840*lambda*cp6*u + 1280*lambda*cp7*u + 10240*sn*m*omega*cp5*up2 - 2560*sn*m*omega*cp6*up2)
    coeffs[40] = ComplexF64(-1152*cp5*u + 2304*cp6*u + 384*cp7*u - 12288*o2*cp4*up3 + 9216*o2*cp5*up3 - 768*m2*cp4*up2 - 384*m2*cp5*up2 + 384*m2*cp6*up2 + 384*lambda*cp4*up2 - 1536*lambda*cp5*up2 + 6912*sn*m*omega*cp4*up2 - 4608*sn*m*omega*cp5*up2)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_6!(coeffs, lambda, m, omega, m2, o2, sn, cp2, cp3, cp4, cp5, cp6, cp7, cp8, u, up2, up3, up4)
    coeffs[41] = ComplexF64(448*cp5*up2 - 3584*o2*cp3*up4 - 448*m2*cp3*up3 - 448*m2*cp4*up3 - 448*lambda*cp4*up2 + 2688*sn*m*omega*cp3*up3)
    coeffs[42] = ComplexF64(-512*o2*cp2*up4 - 128*m2*cp2*up3 - 128*m2*cp3*up3 + 512*sn*m*omega*cp2*up3)
    coeffs[43] = coeffs[28]
    coeffs[44] = ComplexF64(-768*cp7 + 30720*o2*cp7*up2 - 12288*o2*cp8*up2 - 1536*lambda*cp7 + 2304*lambda*cp8 - 4608*sn*m*omega*cp7*u)
    coeffs[45] = ComplexF64(-2688*cp6 + 3840*cp7 + 384*cp8 + 52224*o2*cp6*up2 - 49152*o2*cp7*up2 + 6144*o2*cp8*up2 + 384*m2*cp6*u + 384*m2*cp7*u - 3072*lambda*cp6*u + 4608*lambda*cp7*u - 13824*sn*m*omega*cp6*u + 9216*sn*m*omega*cp7*u)
    coeffs[46] = ComplexF64(-1280*cp5*u + 7680*cp6*u + 1280*cp7*u + 51200*o2*cp5*up3 - 20480*o2*cp6*up3 + 1280*m2*cp5*up2 + 1280*m2*cp6*up2 - 2560*lambda*cp5*u + 7680*lambda*cp6*u - 2560*lambda*cp7*u - 20480*sn*m*omega*cp5*up2 + 5120*sn*m*omega*cp6*up2)
    coeffs[47] = ComplexF64(2880*cp5*u - 5760*cp6*u - 960*cp7*u + 30720*o2*cp4*up3 - 23040*o2*cp5*up3 + 1920*m2*cp4*up2 + 960*m2*cp5*up2 - 960*m2*cp6*up2 - 960*lambda*cp4*up2 + 3840*lambda*cp5*up2 - 17280*sn*m*omega*cp4*up2 + 11520*sn*m*omega*cp5*up2)
    coeffs[48] = ComplexF64(-1344*cp5*up2 + 10752*o2*cp3*up4 + 1344*m2*cp3*up3 + 1344*m2*cp4*up3 + 1344*lambda*cp4*up2 - 8064*sn*m*omega*cp3*up3)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_7!(coeffs, lambda, m, omega, m2, o2, sn, cp2, cp3, cp4, cp5, cp6, cp7, cp8, u, up2, up3, up4)
    coeffs[49] = ComplexF64(1792*o2*cp2*up4 + 448*m2*cp2*up3 + 448*m2*cp3*up3 - 1792*sn*m*omega*cp2*up3)
    coeffs[50] = ComplexF64(256*cp7 - 10240*o2*cp7*up2 + 4096*o2*cp8*up2 + 512*lambda*cp7 - 768*lambda*cp8 + 1536*sn*m*omega*cp7*u)
    coeffs[51] = coeffs[38]
    coeffs[52] = ComplexF64(1280*cp5*u - 7680*cp6*u - 1280*cp7*u - 51200*o2*cp5*up3 + 20480*o2*cp6*up3 - 1280*m2*cp5*up2 - 1280*m2*cp6*up2 + 2560*lambda*cp5*u - 7680*lambda*cp6*u + 2560*lambda*cp7*u + 20480*sn*m*omega*cp5*up2 - 5120*sn*m*omega*cp6*up2)
    coeffs[53] = ComplexF64(-3840*cp5*u + 7680*cp6*u + 1280*cp7*u - 40960*o2*cp4*up3 + 30720*o2*cp5*up3 - 2560*m2*cp4*up2 - 1280*m2*cp5*up2 + 1280*m2*cp6*up2 + 1280*lambda*cp4*up2 - 5120*lambda*cp5*up2 + 23040*sn*m*omega*cp4*up2 - 15360*sn*m*omega*cp5*up2)
    coeffs[54] = ComplexF64(2240*cp5*up2 - 17920*o2*cp3*up4 - 2240*m2*cp3*up3 - 2240*m2*cp4*up3 - 2240*lambda*cp4*up2 + 13440*sn*m*omega*cp3*up3)
    coeffs[55] = ComplexF64(-3584*o2*cp2*up4 - 896*m2*cp2*up3 - 896*m2*cp3*up3 + 3584*sn*m*omega*cp2*up3)
    coeffs[56] = coeffs[30]
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_8!(coeffs, lambda, m, omega, m2, o2, sn, cp2, cp3, cp4, cp5, cp6, cp7, u, up2, up3, up4)
    coeffs[57] = ComplexF64(-640*cp5*u + 3840*cp6*u + 640*cp7*u + 25600*o2*cp5*up3 - 10240*o2*cp6*up3 + 640*m2*cp5*up2 + 640*m2*cp6*up2 - 1280*lambda*cp5*u + 3840*lambda*cp6*u - 1280*lambda*cp7*u - 10240*sn*m*omega*cp5*up2 + 2560*sn*m*omega*cp6*up2)
    coeffs[58] = coeffs[47]
    coeffs[59] = ComplexF64(-2240*cp5*up2 + 17920*o2*cp3*up4 + 2240*m2*cp3*up3 + 2240*m2*cp4*up3 + 2240*lambda*cp4*up2 - 13440*sn*m*omega*cp3*up3)
    coeffs[60] = ComplexF64(4480*o2*cp2*up4 + 1120*m2*cp2*up3 + 1120*m2*cp3*up3 - 4480*sn*m*omega*cp2*up3)
    coeffs[61] = ComplexF64(128*cp5*u - 768*cp6*u - 128*cp7*u - 5120*o2*cp5*up3 + 2048*o2*cp6*up3 - 128*m2*cp5*up2 - 128*m2*cp6*up2 + 256*lambda*cp5*u - 768*lambda*cp6*u + 256*lambda*cp7*u + 2048*sn*m*omega*cp5*up2 - 512*sn*m*omega*cp6*up2)
    coeffs[62] = coeffs[40]
    coeffs[63] = ComplexF64(1344*cp5*up2 - 10752*o2*cp3*up4 - 1344*m2*cp3*up3 - 1344*m2*cp4*up3 - 1344*lambda*cp4*up2 + 8064*sn*m*omega*cp3*up3)
    coeffs[64] = coeffs[55]
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_9!(coeffs, lambda, m, omega, m2, o2, sn, cp3, cp4, cp5, up2, up3, up4)
    coeffs[65] = coeffs[32]
    coeffs[66] = ComplexF64(-448*cp5*up2 + 3584*o2*cp3*up4 + 448*m2*cp3*up3 + 448*m2*cp4*up3 + 448*lambda*cp4*up2 - 2688*sn*m*omega*cp3*up3)
    coeffs[67] = coeffs[49]
    coeffs[68] = ComplexF64(64*cp5*up2 - 512*o2*cp3*up4 - 64*m2*cp3*up3 - 64*m2*cp4*up3 - 64*lambda*cp4*up2 + 384*sn*m*omega*cp3*up3)
    coeffs[69] = coeffs[42]
    coeffs[70] = coeffs[34]
    coeffs[71] = coeffs[11]
    coeffs[72] = coeffs[11]
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_10!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    coeffs[73] = coeffs[11]
    coeffs[74] = coeffs[11]
    coeffs[75] = ComplexF64(256*cp8)
    coeffs[76] = ComplexF64(512*cp7 - 1024*cp8)
    coeffs[77] = ComplexF64(512*cp6 - 1792*cp7 + 1536*cp8)
    coeffs[78] = ComplexF64(256*cp5*up2 - 1024*cp6*up2)
    coeffs[79] = ComplexF64(64*cp4*u - 576*cp5*u + 1024*cp6*u - 256*cp7*u)
    coeffs[80] = ComplexF64(-128*cp4*up2 + 256*cp5*up2)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_11!(coeffs, cp4, cp5, cp6, cp7, cp8, up2)
    coeffs[81] = ComplexF64(64*cp4*up2)
    coeffs[82] = coeffs[11]
    coeffs[83] = coeffs[11]
    coeffs[84] = coeffs[11]
    coeffs[85] = ComplexF64(-1024*cp8)
    coeffs[86] = ComplexF64(-2560*cp7 + 5120*cp8)
    coeffs[87] = ComplexF64(-3072*cp6 + 10752*cp7 - 9216*cp8)
    coeffs[88] = ComplexF64(-1792*cp5*up2 + 7168*cp6*up2)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_12!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    coeffs[89] = ComplexF64(-512*cp4*u + 4608*cp5*u - 8192*cp6*u + 2048*cp7*u)
    coeffs[90] = ComplexF64(1152*cp4*up2 - 2304*cp5*up2)
    coeffs[91] = ComplexF64(-640*cp4*up2)
    coeffs[92] = coeffs[11]
    coeffs[93] = coeffs[11]
    coeffs[94] = ComplexF64(1536*cp8)
    coeffs[95] = ComplexF64(5120*cp7 - 10240*cp8)
    coeffs[96] = ComplexF64(7680*cp6 - 26880*cp7 + 23040*cp8)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_13!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    coeffs[97] = ComplexF64(5376*cp5*up2 - 21504*cp6*up2)
    coeffs[98] = ComplexF64(1792*cp4*u - 16128*cp5*u + 28672*cp6*u - 7168*cp7*u)
    coeffs[99] = ComplexF64(-4608*cp4*up2 + 9216*cp5*up2)
    coeffs[100] = ComplexF64(2880*cp4*up2)
    coeffs[101] = coeffs[11]
    coeffs[102] = coeffs[85]
    coeffs[103] = ComplexF64(-5120*cp7 + 10240*cp8)
    coeffs[104] = ComplexF64(-10240*cp6 + 35840*cp7 - 30720*cp8)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_14!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    coeffs[105] = ComplexF64(-8960*cp5*up2 + 35840*cp6*up2)
    coeffs[106] = ComplexF64(-3584*cp4*u + 32256*cp5*u - 57344*cp6*u + 14336*cp7*u)
    coeffs[107] = ComplexF64(10752*cp4*up2 - 21504*cp5*up2)
    coeffs[108] = ComplexF64(-7680*cp4*up2)
    coeffs[109] = coeffs[75]
    coeffs[110] = ComplexF64(2560*cp7 - 5120*cp8)
    coeffs[111] = coeffs[96]
    coeffs[112] = ComplexF64(8960*cp5*up2 - 35840*cp6*up2)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_15!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    coeffs[113] = ComplexF64(4480*cp4*u - 40320*cp5*u + 71680*cp6*u - 17920*cp7*u)
    coeffs[114] = ComplexF64(-16128*cp4*up2 + 32256*cp5*up2)
    coeffs[115] = ComplexF64(13440*cp4*up2)
    coeffs[116] = ComplexF64(-512*cp7 + 1024*cp8)
    coeffs[117] = coeffs[87]
    coeffs[118] = ComplexF64(-5376*cp5*up2 + 21504*cp6*up2)
    coeffs[119] = coeffs[106]
    coeffs[120] = ComplexF64(16128*cp4*up2 - 32256*cp5*up2)
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_16!(coeffs, cp4, cp5, cp6, up2)
    coeffs[121] = ComplexF64(-16128*cp4*up2)
    coeffs[122] = coeffs[77]
    coeffs[123] = ComplexF64(1792*cp5*up2 - 7168*cp6*up2)
    coeffs[124] = coeffs[98]
    coeffs[125] = ComplexF64(-10752*cp4*up2 + 21504*cp5*up2)
    coeffs[126] = coeffs[115]
    coeffs[127] = ComplexF64(-256*cp5*up2 + 1024*cp6*up2)
    coeffs[128] = coeffs[89]
    return nothing
end

@noinline function _construct_gsn_basis_yc_s0_block_17!(coeffs, cp4, cp5, up2)
    coeffs[129] = ComplexF64(4608*cp4*up2 - 9216*cp5*up2)
    coeffs[130] = coeffs[108]
    coeffs[131] = coeffs[79]
    coeffs[132] = ComplexF64(-1152*cp4*up2 + 2304*cp5*up2)
    coeffs[133] = coeffs[100]
    coeffs[134] = ComplexF64(128*cp4*up2 - 256*cp5*up2)
    coeffs[135] = coeffs[91]
    coeffs[136] = coeffs[81]
    return nothing
end

function construct_gsn_basis_yc_s0(lambda, m, nu, omega)
    lambda = ComplexF64(lambda)
    m = Float64(m)
    nu = Float64(nu)
    omega = ComplexF64(omega)
    m2 = m * m
    o2 = omega * omega
    s1 = sin(nu)
    c1 = sqrt((1 - s1) * (1 + s1))
    sn = s1
    c = c1
    cp2 = c * c
    cp3 = cp2 * c
    cp4 = cp3 * c
    cp5 = cp4 * c
    cp6 = cp5 * c
    cp7 = cp6 * c
    cp8 = cp7 * c
    cp9 = cp8 * c
    cp10 = cp9 * c
    u = s1 * s1 / (1 + c1)
    up2 = u * u
    up3 = up2 * u
    up4 = up3 * u
    coeffs = Vector{ComplexF64}(undef, LEN_GSN_BASIS_yc_s0)
    _construct_gsn_basis_yc_s0_block_1!(coeffs, cp2, cp3, cp4, u)
    _construct_gsn_basis_yc_s0_block_2!(coeffs, cp2, cp3, cp4, u)
    _construct_gsn_basis_yc_s0_block_3!(coeffs, cp2, cp3, cp4, u)
    _construct_gsn_basis_yc_s0_block_4!(coeffs, lambda, m, omega, m2, o2, sn, cp4, cp5, cp6, cp7, cp8, cp9, cp10, u, up2, up3)
    _construct_gsn_basis_yc_s0_block_5!(coeffs, lambda, m, omega, m2, o2, sn, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9, u, up2, up3, up4)
    _construct_gsn_basis_yc_s0_block_6!(coeffs, lambda, m, omega, m2, o2, sn, cp2, cp3, cp4, cp5, cp6, cp7, cp8, u, up2, up3, up4)
    _construct_gsn_basis_yc_s0_block_7!(coeffs, lambda, m, omega, m2, o2, sn, cp2, cp3, cp4, cp5, cp6, cp7, cp8, u, up2, up3, up4)
    _construct_gsn_basis_yc_s0_block_8!(coeffs, lambda, m, omega, m2, o2, sn, cp2, cp3, cp4, cp5, cp6, cp7, u, up2, up3, up4)
    _construct_gsn_basis_yc_s0_block_9!(coeffs, lambda, m, omega, m2, o2, sn, cp3, cp4, cp5, up2, up3, up4)
    _construct_gsn_basis_yc_s0_block_10!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    _construct_gsn_basis_yc_s0_block_11!(coeffs, cp4, cp5, cp6, cp7, cp8, up2)
    _construct_gsn_basis_yc_s0_block_12!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    _construct_gsn_basis_yc_s0_block_13!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    _construct_gsn_basis_yc_s0_block_14!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    _construct_gsn_basis_yc_s0_block_15!(coeffs, cp4, cp5, cp6, cp7, cp8, u, up2)
    _construct_gsn_basis_yc_s0_block_16!(coeffs, cp4, cp5, cp6, up2)
    _construct_gsn_basis_yc_s0_block_17!(coeffs, cp4, cp5, up2)
    return coeffs
end

function eval_gsn_basis_yc_s0!(out::Vector{ComplexF64}, coeffs::Vector{ComplexF64}, t::ComplexF64)
    return _evaluate_polynomial_blocks!(out,coeffs,t,(1:4,5:7,8:9,10:10,1:0,1:0,1:0,1:0,1:0,1:0,11:15,16:19,20:22,23:24,25:25,1:0,1:0,1:0,1:0,1:0,1:0,26:34,35:42,43:49,50:55,56:60,61:64,65:67,68:69,70:70,71:81,82:91,92:100,101:108,109:115,116:121,122:126,127:130,131:133,134:135,136:136,1:0,1:0,))
end

ordinary_value_count() = 43

function construct_gsn_basis(lambda, m, nu, omega)
    return construct_gsn_basis_yc_s0(lambda, m, nu, omega)
end

function eval_gsn_basis!(out::Vector{ComplexF64}, coeffs::Vector{ComplexF64}, x0)
    x0 = ComplexF64(x0)
    return eval_gsn_basis_yc_s0!(out, coeffs, 1.0 - x0)
end

function eval_gsn_basis_variable!(out::Vector{ComplexF64}, coeffs::Vector{ComplexF64}, value)
    return eval_gsn_basis_yc_s0!(out, coeffs, ComplexF64(value))
end

end
