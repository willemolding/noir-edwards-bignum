use crate::Curve;
use crate::TECurveParameterTrait;
use bignum::BigNum;
use bignum::ED25519_Fq;

/// Parameters for the ed25519 twisted Edwards curve.
///
/// ed25519 is defined over the prime field `2^255 - 19` (provided by `ED25519_Fq`) with:
///   a = -1
///   d = -121665/121666 mod p
/// and the standard base point B with `y = 4/5`.
///
pub struct Ed25519Params {}
impl TECurveParameterTrait<ED25519_Fq> for Ed25519Params {
    fn a() -> ED25519_Fq {
        // a = -1 mod p
        BigNum::from_limbs([
            0xffffffffffffffffffffffffffffec,
            0xffffffffffffffffffffffffffffff,
            0x7fff,
        ])
    }
    fn d() -> ED25519_Fq {
        // d = -121665/121666 mod p
        BigNum::from_limbs([
            0x700a4d4141d8ab75eb4dca135978a3,
            0x6cee2b6ffe738cc740797779e89800,
            0x5203,
        ])
    }
    fn gen() -> (ED25519_Fq, ED25519_Fq) {
        (
            BigNum::from_limbs([
                0x2cc7609525a7b2c9562d608f25d51a,
                0x36d3cd6e53fec0a4e231fdd6dc5c69,
                0x2169,
            ]),
            BigNum::from_limbs([
                0x666666666666666666666666666658,
                0x666666666666666666666666666666,
                0x6666,
            ]),
        )
    }
}
pub type Ed25519 = Curve<ED25519_Fq, Ed25519Params>;
