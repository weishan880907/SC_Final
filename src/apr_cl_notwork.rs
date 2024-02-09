// 501900 Ximing Zhang

use num_bigint::BigInt;
use num_traits::ToPrimitive;
use num_integer::Integer;
use std::collections::HashMap;
use std::fmt;
use std::collections::HashSet;
use std::ops::Mul;
use std::{env, fs};
use std::str::FromStr;
use std::time::Instant;

fn isprime_slow(n: &mut BigInt) -> bool {
    if n < &mut BigInt::from(2) {
        return false;

    } else if  n == &mut BigInt::from(2) ||  n == &mut BigInt::from(3) {
        return true;

    } else if  n.clone() % BigInt::from(2) == BigInt::from(0) {
        return false;

    } else {

        let mut i = BigInt::from(3);
        while &mut (i.clone() * i.clone()) <= n {

            if n.clone() % i.clone() == BigInt::from(0) {

                return false;
            }

            i += BigInt::from(2);
        }
    }
    return true;
}

fn v(qq: &mut BigInt, tt: &mut BigInt) -> u32 {
    let mut ans = BigInt::from(0);
    let  q = qq.clone();
    let mut t = tt.clone();


    while t.clone() % q.clone() == BigInt::from(0) {
        ans += BigInt::from(1);
        t /= &q;
    }
    let mut res:u32 = 0;

    match ans.to_u32() {
        Some(i) => {
            res = i;
        },
        None => {
            // println!("Unable to convert BigInt to i64");
        }
    }
    return res;
}

fn e(tt: &mut BigInt) -> (BigInt, Vec<BigInt>) {
    let mut t = tt.clone();
    let mut s = BigInt::from(1);
    let mut q_list:Vec<BigInt> = Vec::new();
    let  start:i64 = 2;
    let mut end:i64 = 2;


    //BigInt convert to i64
    match t.to_i64() {
        Some(i) => {
            end = i + 2;
            // println!("end: {}", end);
        },
        None => {
            // println!("Unable to convert BigInt to i64");
        }
    }


    // println!("{},{}", start, end);
    for qt in start..end {
        let mut q = BigInt::from(qt);
        if t.clone() % (q.clone() - BigInt::from(1)) == BigInt::from(0) && isprime_slow(&mut q){
            // println!("{},{}", q, t);

            let  mi:u32 = 1 + v(&mut q,&mut t);
            // println!("mi:{}", mi);
            s *= q.pow(mi);
            q_list.push(q.clone());
        }
    }

    let  res = BigInt::from(2) * s;
    let  res_q_list = q_list.clone();
    return (res.clone(), res_q_list.clone());
}

fn prime_factorize(n: &mut BigInt) -> Vec<(BigInt,BigInt)> {
    let mut ret:Vec<(BigInt,BigInt)> = Vec::new();
    let mut p = BigInt::from(2);

    while &mut (p.clone() * p.clone()) <= n {

        if n.clone() % p.clone() == BigInt::from(0) {

            let mut num = BigInt::from(0);
            while (n.clone() % p.clone()) == BigInt::from(0) {

                num += BigInt::from(1);
                *n /= p.clone();
            }
            ret.push((p.clone(), num.clone()));
        }
        p += BigInt::from(1);
    }
    if n.clone() != BigInt::from(1) {
        ret.push((n.clone(), BigInt::from(1)));
    }
    return ret;
}
#[derive(Clone)]
struct JacobiSum {
    p: BigInt,
    k: BigInt,
    q: BigInt,
    m: BigInt,
    pk: BigInt,
    coef: Vec<BigInt>,
}
// Jacob element multiply Jacob element
impl Mul<JacobiSum> for JacobiSum {
    type Output = Self;
    fn mul(mut self, jac: JacobiSum) -> Self::Output {

        let  m = self.m.clone();
        let  pk = self.pk.clone();
        let mut j_ret = JacobiSum::new(&mut self.p,&mut self.k,&mut self.q);
        // println!("{:?},{:?},{:?}",m,pk,j_ret);
        let mut temp_m:u32 = 0;
        match m.to_u32() {
            Some(i) => {
                temp_m = i;
            },
            None => {
                // println!("temp_m is unable to convert BigInt to i64");
            }
        }
        for i in 0..temp_m {

            for j in 0..temp_m {

                let  ii = BigInt::from(i);
                let  jj = BigInt::from(j);
                if (ii.clone() + jj.clone()) % pk.clone() < m {
                    let  temp = (ii.clone() + jj.clone()) % pk.clone();
                    let mut temp_ii:usize = 0;
                    match temp.to_usize() {
                        Some(t) => {
                            temp_ii = t;
                        },
                        None => {
                            // println!("temp_m is unable to convert BigInt to i64");
                        }
                    }
                    // let mut ui:usize = i as usize;
                    j_ret.coef[temp_ii] += self.coef[i as usize].clone() * jac.coef[j as usize].clone();

                } else {
                    let mut temp_k:u32 = 0;
                    match self.k.to_u32() {
                        Some(t) => {
                            temp_k = t;
                        },
                        None => {
                            // println!("temp_m is unable to convert BigInt to i64");
                        }
                    }
                    let mut r = (ii.clone() + jj.clone()) % pk.clone() - self.p.pow(temp_k - 1);
                    while r >= BigInt::from(0) {
                        let mut temp_r:usize = 0;
                        match r.to_usize() {
                            Some(i) => {
                                temp_r = i;
                            },
                            None => {
                                // println!("temp_m is unable to convert BigInt to i64");
                            }
                        }
                        j_ret.coef[temp_r] -= self.coef[i as usize].clone() * jac.coef[j as usize].clone(); //encounter negative number 
                        // j_ret.coef[temp_r] = j_ret.coef[temp_r].clone() + 
                        // println!("{:?}", j_ret.coef[temp_r]);
                        r -= self.p.pow(temp_k - 1);
                    }
                }
            }
        }
        return j_ret.clone();
    }
}

//Jacob element mutiply  BigInt
impl Mul<BigInt> for JacobiSum {
    
    type Output = Self;
    fn mul(mut self, rhs: BigInt) -> Self::Output {

        let  right = rhs.clone();
        let mut j_ret = JacobiSum::new(&mut self.p,&mut self.k,&mut self.q);
        let mut temp_m:u32 = 0;
        match self.m.to_u32() {
            Some(i) => {
                temp_m = i;
            },
            None => {
                println!("temp_m is unable to convert BigInt to i64");
            }
        }
        for i in 0..temp_m {
            j_ret.coef[i as usize] = self.coef[i as usize].clone() * right.clone();
        }
        return j_ret.clone();
        // unimplemented!()
    }
}


impl fmt::Debug for JacobiSum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("JacobiSum")
            .field("p", &self.p)
            .field("k", &self.k)
            .field("q", &self.q)
            .field("m", &self.m)
            .field("pk", &self.pk)
            .field("coef", &self.coef)
            .finish()
    }
}


impl JacobiSum {
    
    fn new(pp: &mut BigInt, kk: &mut BigInt, qq: &mut BigInt) -> Self {

        let  pt = pp.clone();
        let  kt = kk.clone();
        let  qt = qq.clone();
        let mut temp_k:u32 = 0;
        match kt.to_u32() {
            Some(i) => {
                temp_k = i;
            },
            None => {
                // println!("temp_m is unable to convert BigInt to i64");
            }
        }
        let  mt = (pp.clone() - BigInt::from(1)) * pp.pow(temp_k - 1);
        let  pkt = pp.pow(temp_k);
        let mut temp_m:u32 = 0;
        match mt.to_u32() {
            Some(i) => {
                temp_m = i;
            },
            None => {
                // println!("temp_m is unable to convert BigInt to i64");
            }
        }
        let mut coeft:Vec<BigInt> = Vec::new();
        for _ in 0..temp_m {
            coeft.push(BigInt::from(0));
        }
        JacobiSum{p:pt.clone(), k:kt.clone(), q:qt.clone(), m:mt.clone(), pk:pkt.clone(), coef: coeft.clone()}
    }
    fn mmod(&mut self, n:&mut BigInt) -> Self{

        let  nn = n.clone();
        let mut temp_m:usize = 0;
        match self.m.to_usize() {
            Some(i) => {
                temp_m = i;
            },
            None => {
                // println!("temp_m is unable to convert BigInt to i64");
            }
        }
        for i in 0..temp_m {
            self.coef[i] = (self.coef[i].clone() +  nn.clone()) % nn.clone();
        }
        return self.clone();
    }
    fn modpow(&mut self, x:&mut BigInt, n:&mut BigInt) -> JacobiSum {
        let mut j_ret = JacobiSum::new(&mut self.p,&mut self.k,&mut self.q);
        j_ret.coef[0] = BigInt::from(1);
        let mut j_a = self.clone();
        // println!("{:?}", j_a);
        let mut xx = x.clone();
        let mut nn = n.clone();
        while xx.clone() > BigInt::from(0) {
            if xx.clone() % BigInt::from(2) == BigInt::from(1) {
                j_ret = (j_ret.clone() * j_a.clone()).mmod(&mut nn);
            }
            j_a = j_a.clone() * j_a.clone();
            j_a = j_a.mmod(&mut nn);
            // println!("{:?}", j_a);
            xx /= BigInt::from(2);
        }
        // println!("{:?}", j_ret);
        j_ret = j_ret.mmod(&mut nn);
        return j_ret.clone();
    }
    fn one(&mut self) -> Self { 
        self.coef[0] = BigInt::from(1);
        let mut temp_m:usize = 0;
        match self.m.to_usize() {
            Some(i) => {
                temp_m = i;
            },
            None => {
                println!("temp_m is unable to convert BigInt to i64");
            }
        }
        for i in 1..temp_m {
            self.coef[i] = BigInt::from(0);
        }
        return self.clone();
    }
    fn sigma_inv(&mut self, x: &mut BigInt) -> JacobiSum {
        let  m = self.m.clone();
        let  pk = self.pk.clone();
        let mut j_ret = JacobiSum::new(&mut self.p,&mut self.k,&mut self.q);
        let mut temp_pk:u32 = 0;
        match pk.to_u32() {
            Some(i) => {
                temp_pk = i;
            },
            None => {
                println!("temp_m is unable to convert BigInt to i64");
            }
        }
        // let mut xx = x.clone();
        for i in 0..temp_pk {
            let  ii = BigInt::from(i);
            if ii < m {
                if (ii.clone() * x.clone()) % pk.clone() < m {
                    let mut temp_ii:usize = 0;
                    match ii.to_usize() {
                        Some(i) => {
                            temp_ii = i;
                        },
                        None => {
                            println!("temp_m is unable to convert BigInt to i64");
                        }
                    }
                    let  jj = (ii.clone() * x.clone()) % pk.clone();
                    let mut temp_jj:usize = 0;
                    match jj.to_usize() {
                        Some(i) => {
                            temp_jj = i;
                        },
                        None => {
                            // println!("temp_m is unable to convert BigInt to i64");
                        }
                    }
                    j_ret.coef[temp_ii] += self.coef[temp_jj].clone();
                }
            } else {
                let mut temp_k:u32 = 0;
                match self.k.to_u32() {
                    Some(i) => {
                        temp_k = i;
                    },
                    None => {
                        // println!("temp_m is unable to convert BigInt to i64");
                    }
                }
                let mut r = ii.clone() - self.p.pow(temp_k - 1);
                while r >= BigInt::from(0) {
                    if (ii.clone() * x.clone()) % pk.clone() < m {
                        let mut temp_ii:usize = 0;
                        match r.to_usize() {
                            Some(i) => {
                                temp_ii = i;
                            },
                            None => {
                                // println!("temp_m is unable to convert BigInt to i64");
                            }
                        }
                        let  jj = (ii.clone() * x.clone()) % pk.clone();
                        let mut temp_jj:usize = 0;
                        match jj.to_usize() {
                            Some(i) => {
                                temp_jj = i;
                            },
                            None => {
                                // println!("temp_m is unable to convert BigInt to i64");
                            }
                        }
                        j_ret.coef[temp_ii] -= self.coef[temp_jj].clone();
                    }
                    r -= self.p.pow(temp_k - 1);
                }
            }
        }
        return j_ret.clone();
    }

    fn is_root_of_unity(&mut self, n: &mut BigInt) -> (bool, BigInt) {
        let  n = n.clone();
        let  m = self.m.clone();
        let  p = self.p.clone();
        let  k = self.k.clone();
        let mut one = BigInt::from(0);
        let mut temp_m:usize = 0;
        match m.to_usize() {
            Some(i) => {
                temp_m = i;
            },
            None => {
                // println!("temp_m无法将 BigInt 转换为 i64");
            }
        }
        let mut h = BigInt::from(0);
        for i in 0..temp_m {
            if self.coef[i] == BigInt::from(1) {
                one += BigInt::from(1);
                h = BigInt::from(i);
            } else if self.coef[i] == BigInt::from(0) {
                continue;
            } else if (self.coef[i].clone() - BigInt::from(-1)) % n.clone() != BigInt::from(0) {
                return (false, BigInt::from(0).clone());// BigInt::from(-1).clone()表示None
            }
        }
        if one == BigInt::from(1) {
            return (true, h.clone());
        }
        let mut temp_i = BigInt::from(0);
        for i in 0..temp_m {
            if self.coef[i] != BigInt::from(0) {
                temp_i = BigInt::from(i);
                break;
            }
        }
        let mut temp_k:u32 = 0;
        match k.to_u32() {
            Some(i) => {
                temp_k = i;
            },
            None => {
                // println!("temp_m is unable to convert BigInt to i64");
            }
        }
        let  r = temp_i.clone() % (p.pow(temp_k - 1));
        for i in 0..temp_m {
            let  ii = BigInt::from(i);
            if ii.clone() % (p.pow(temp_k - 1)) == r {
                if (self.coef[i].clone() - BigInt::from(-1)) % n.clone() != BigInt::from(0) {
                    return (false, BigInt::from(0).clone());// BigInt::from(-1).clone()表示None
                }
            } else {
                if self.coef[i] != BigInt::from(0) {
                    return (false, BigInt::from(0).clone());// BigInt::from(-1).clone()表示None 
                }
            }
        }
        let  temp = (p.clone() - BigInt::from(1)) * p.pow(temp_k - 1) + r.clone();
        return (true, temp.clone());
    }

}

fn smallest_primitive_root(q: &mut BigInt) -> BigInt {
    let  qq = q.clone();
    let mut temp_q:u32 = 0;
    match qq.to_u32() {
        Some(i) => {
            temp_q = i;
        },
        None => {
            // println!("smallest_primitive_root is unable to convert BigInt to i64");
        }
    }
    for r in 2..temp_q {
        let mut s: HashSet<BigInt> = HashSet::new();
        let mut m = BigInt::from(1);
        let  rr = BigInt::from(r);
        for _ in 1..temp_q {
            m = (m.clone() * rr.clone()) % qq.clone();
            s.insert(m.clone());
        }
        let  set_len = BigInt::from(s.len());
        if set_len == qq.clone() - BigInt::from(1) {
            return rr.clone();
        }
    }
    return BigInt::from(0).clone(); 
}


fn calc_f(q: &mut BigInt) -> HashMap<BigInt, BigInt> {
    let mut qq = q.clone();
    let  g = smallest_primitive_root(&mut qq);
    let mut m: HashMap<BigInt, BigInt> = HashMap::new();
    let mut temp_q:u32 = 0;
    match qq.to_u32() {
        Some(i) => {
            temp_q = i;
        },
        None => {
            // println!("calc_f is unable to convert BigInt to i64");
        }
    }
    for x in 1..temp_q - 1 {
        let  k = g.pow(x) % q.clone();
        let  xx = BigInt::from(x);
        m.insert(k.clone(), xx.clone());
    }
    let mut f: HashMap<BigInt, BigInt> = HashMap::new();
    for x in 1..temp_q - 1 {
        
        let  m_k = (BigInt::from(1) - (g.pow(x) % q.clone()) + q.clone()) % q.clone();
        // println!("m_k {:?}", m_k);
        let  f_k = BigInt::from(x);
        let key = &m_k;
        let mut fk = BigInt::from(0);
        match m.get(key) {
            Some(value) => {
                fk = value.clone();
            }
            None => {
                // println!("No value found for key: {}", key)
            },
        }
        f.insert(f_k.clone(), fk.clone());
    }
    return f.clone();
}





fn calc_j_ab(p: &mut BigInt,k: &mut BigInt,q: &mut BigInt, a: &mut BigInt, b: &mut BigInt) -> JacobiSum {
    let mut pp = p.clone();
    let mut kk = k.clone();
    let mut qq = q.clone();
    let  aa = a.clone();
    let  bb = b.clone();
    let mut j_ret = JacobiSum::new(&mut pp,&mut kk, &mut qq);
    let  f = calc_f(&mut qq);
    // println!("{:?}", f);
    let mut temp_q:u32 = 0;
    match qq.to_u32() {
        Some(i) => {
            temp_q = i;
        },
        None => {
            // println!("calc_J_ab is unable to convert BigInt to i64");
        }
    }
    let mut temp_k:u32 = 0;
    match kk.to_u32() {
        Some(i) => {
            temp_k = i;
        },
        None => {
            // println!("calc_J_ab is unable to convert BigInt to i64");
        }
    }
    for x in 1..temp_q - 1{
        let  pk = pp.pow(temp_k);
        let  xx = BigInt::from(x);
        let mut k_xx = BigInt::from(0);
        match f.get(&xx) {
            Some(value) => {
                k_xx = value.clone();
            }
            None => {
                println!("Key {} not found", xx);
            }
        }
        let  temp = (aa.clone() * xx.clone() + bb.clone() * k_xx) % pk.clone();
        if temp < j_ret.m {
            let mut temp_temp:usize = 0;
            match temp.to_usize() {
                Some(i) => {
                    temp_temp = i;
                },
                None => {
                    // println!("calc_J_ab temp_temp is unable to convert BigInt to i64");
                }
            }
            j_ret.coef[temp_temp] += BigInt::from(1);
        } else {
            let mut r = temp - pp.pow(temp_k - 1);
            while r >= BigInt::from(0) {
                let mut temp_r:usize = 0;
                match r.to_usize() {
                    Some(i) => {
                        temp_r = i;
                    },
                    None => {
                        // println!("calc_J_ab temp_r is unable to convert BigInt to i64");
                    }
                }
                j_ret.coef[temp_r] -= BigInt::from(1);
                r -= pp.pow(temp_k - 1);
            }
        }
    }
    return j_ret.clone();
}

fn calc_j(p: &mut BigInt,k: &mut BigInt,q: &mut BigInt) -> JacobiSum{
    let mut pp = p.clone();
    let mut kk = k.clone();
    let mut qq = q.clone();
    return calc_j_ab(&mut pp,&mut kk,&mut qq, &mut BigInt::from(1), &mut BigInt::from(1)).clone();
}


fn aprtest_step4a(p: &mut BigInt,k: &mut BigInt, q: &mut BigInt, n: &mut BigInt) -> (bool, BigInt) {
    let mut pp = p.clone();
    let mut kk = k.clone();
    let mut qq = q.clone();
    let mut nn = n.clone();
    //println!("Step 4a. (p^k, q = {0}^{1}, {2})",pp, kk, qq);
    let mut j = calc_j(&mut pp,&mut kk,&mut qq);
    // println!("{:?})", J);
    let mut s1 = JacobiSum::new(&mut pp,&mut kk, &mut qq).one();
    let mut temp_k:u32 = 0;
    match kk.to_u32() {
        Some(i) => {
            temp_k = i;
        },
        None => {
            // println!("calc_J_ab temp_r  is unable to convert BigInt to i64");
        }
    }
    let  up = pp.pow(temp_k);
    let mut temp_up:u32 = 0;
    match up.to_u32() {
        Some(i) => {
            temp_up = i;
        },
        None => {
            // println!("APRtest_step4b temp_k is unable to convert BigInt to i64");
        }
    }
    // println!("{:?}", temp_up);
    for x in 0..temp_up {
        let mut xx = BigInt::from(x);
        if xx.clone() % pp.clone() == BigInt::from(0) {
            continue;
        }
        let mut t = j.sigma_inv(&mut xx);
        // println!("{:?}, {:?}", xx, NN);
        t = t.modpow(&mut xx, &mut nn);
        // println!("{:?}", t);
        s1 = s1 * t;
        s1 = s1.mmod(&mut nn);
    }
    s1 = s1.mmod(&mut nn);
    let  r = nn.clone() % (pp.pow(temp_k));
    let mut can = nn.clone() / pp.pow(temp_k);
    // println!("{:?}, {:?}, {:?}",s1 , can, NN);
    let  s2 = s1.modpow(&mut can, &mut nn);
    // println!("{:?}", s2);
    let mut j_alpha = JacobiSum::new(&mut pp,&mut kk, &mut qq).one();
    let mut temp_pk:u32 = 0;
    let  temp = pp.pow(temp_k);
    match temp.to_u32() {
        Some(i) => {
            temp_pk = i;
        },
        None => {
            // println!("APRtest_step4a is unable to convert BigInt to i64");
        }
    }
    for x in 0..temp_pk {
        let mut xx = BigInt::from(x);
        if xx.clone() % pp.clone() == BigInt::from(0) {
            continue;
        }
        let mut t = j.sigma_inv(&mut xx);
        t = t.modpow(&mut ((r.clone() * xx.clone())/pp.pow(temp_k)), &mut nn);
        j_alpha = j_alpha.clone() * t.clone();
        j_alpha = j_alpha.mmod(&mut nn);
    }

    let mut s = (s2.clone() * j_alpha.clone()).mmod(&mut nn);
    s = s.mmod(&mut nn);
    // println!("{:?}", S);
    let (exist, h) = s.is_root_of_unity(&mut nn);

    if !exist {
        return (false, BigInt::from(0).clone()); // BigInt::from(-1).clone() represent None 
    } else {
        let mut l_p = BigInt::from(0);
        
        if h.clone() % pp.clone() != l_p {
            l_p = BigInt::from(1);
        } else {
            l_p = BigInt::from(0);
        }
        return (true, l_p.clone());
    }
}


fn calc_j3(p:&mut BigInt,k:&mut BigInt,q:&mut BigInt) -> JacobiSum {
    let mut pp = p.clone();
    let mut kk = k.clone();
    let mut qq = q.clone();
    let  j2q = calc_j(&mut pp,&mut kk,&mut qq);
    let  j21 = calc_j_ab(&mut pp,&mut kk,&mut qq, &mut BigInt::from(2), &mut BigInt::from(1));
    let  j_ret = j2q * j21;
    return j_ret.clone();
}

fn calc_j2(p:&mut BigInt,k:&mut BigInt,q:&mut BigInt) -> JacobiSum {
    let mut pp = p.clone();
    let mut kk = k.clone();
    let mut qq = q.clone();
    let j31 = calc_j_ab(&mut BigInt::from(2),&mut BigInt::from(3),&mut qq, &mut BigInt::from(3), &mut BigInt::from(1));
    let mut j_conv = JacobiSum::new(&mut pp,&mut kk, &mut qq);
    let mut temp_m:u32 = 0;
    match j31.m.to_u32() {
        Some(i) => {
            temp_m = i;
        },
        None => {
            // println!("calc_J2 temp_m is unable to convert BigInt to i64");
        }
    }
    let mut temp_k:u32 = 0;
    match kk.to_u32() {
        Some(i) => {
            temp_k = i;
        },
        None => {
            //  println!("calc_J2 temp_k is unable to convert BigInt to i64");
        }
    }
    for i in 0..temp_m {
        let index = BigInt::from(i) * pp.pow(temp_k) / BigInt::from(8);
        let mut temp_index:usize = 0;
        match index.to_usize() {
            Some(i) => {
                temp_index = i;
            },
            None => {
                //  println!("calc_J2 temp_index is unable to convert BigInt to i64");
            }
        }
        j_conv.coef[temp_index] = j31.coef[i as usize].clone();
    }
    let  j_ret = j_conv.clone() * j_conv.clone();
    return j_ret.clone();
}
fn aprtest_step4b(p: &mut BigInt,k: &mut BigInt, q: &mut BigInt, n: &mut BigInt) -> (bool, BigInt) {
    let mut pp = p.clone();
    let mut kk = k.clone();
    let mut qq = q.clone();
    let mut nn = n.clone();
    //println!("Step 4b. (p^k, q = {0}^{1}, {2})",pp, kk, qq);
    let mut j = calc_j3(&mut pp,&mut kk,&mut qq);
    // println!("{:?}", J);
    let mut s1 = JacobiSum::new(&mut pp,&mut kk, &mut qq).one();
    // println!("{:?}", s1);
    let mut temp_k:u32 = 0;
    match kk.to_u32() {
        Some(i) => {
            temp_k = i;
        },
        None => {
            // println!("APRtest_step4b temp_k is unable to convert BigInt to i64");
        }
    }
    let  up = pp.pow(temp_k);
    let mut temp_up:u32 = 0;
    match up.to_u32() {
        Some(i) => {
            temp_up = i;
        },
        None => {
            // println!("APRtest_step4b temp_k is unable to convert BigInt to i64");
        }
    }
    // println!("{:?}", temp_up);
    for x in 0..temp_up {
        let mut xx = BigInt::from(x);
        if xx.clone() % BigInt::from(8) != BigInt::from(1) && xx.clone() % BigInt::from(8) != BigInt::from(3) {
            continue;
        }
        
        let mut t = j.sigma_inv(&mut xx);
        // println!("{:?},{:?}", xx, NN);
        t = t.modpow(&mut xx, &mut nn);
        // println!("{:?}", t);
        s1 = s1 * t;
        s1 = s1.mmod(&mut nn);
    }
    s1 = s1.mmod(&mut nn);
    // println!("{:?}", s1);
    let mut temp_k:u32 = 0;
    match kk.to_u32() {
        Some(i) => {
            temp_k = i;
        },
        None => {
            // println!("APRtest_step4b is unable to convert BigInt to i64");
        }
    }
    let  r = nn.clone() % (pp.pow(temp_k));
    let mut can = nn.clone() / pp.pow(temp_k);
    let  s2 = s1.modpow(&mut can, &mut nn);
    // println!("{:?}", s2);
    let mut j_alpha = JacobiSum::new(&mut pp,&mut kk, &mut qq).one();
    let mut temp_pk:u32 = 0;
    let  temp = pp.pow(temp_k);
    match temp.to_u32() {
        Some(i) => {
            temp_pk = i;
        },
        None => {
            // println!("APRtest_step4b is unable to convert BigInt to i64");
        }
    }
    // println!("{:?}",temp_pk);
    for x in 0..temp_pk {
        let mut xx = BigInt::from(x);
        if xx.clone() % BigInt::from(8) != BigInt::from(1) && xx.clone() % BigInt::from(8) != BigInt::from(3) {
            continue;
        }
        let mut t = j.sigma_inv(&mut xx);
        t = t.modpow(&mut ((r.clone() * xx.clone())/pp.pow(temp_k)), &mut nn);
        j_alpha = j_alpha.clone() * t.clone();
        j_alpha = j_alpha.mmod(&mut nn);
    }
    #[allow(unused_assignments)]
    let mut s = JacobiSum::new(&mut pp,&mut kk, &mut qq);
    // 
    if nn.clone() % BigInt::from(8) == BigInt::from(1) || nn.clone() % BigInt::from(8) == BigInt::from(3) {
        let mut temp = s2.clone() * j_alpha.clone();
        // println!("{:?}", temp);
        s = temp.mmod(&mut nn);
        
    } else {
        let  j2_delta = calc_j2(&mut pp,&mut kk, &mut qq);
        let mut temp = s2.clone() * j_alpha.clone();
        temp = temp.clone() * j2_delta.clone();
        s = temp.mmod(&mut nn);
    }
    s = s.mmod(&mut nn);
    // println!("{:?}, {:?}", S, NN);
    let (exist, h) = s.is_root_of_unity(&mut nn);
    // println!("{:?}, {:?}", exist, h);
    if !exist {
        return (false, BigInt::from(0).clone()); // BigInt::from(-1).clone() represent None 
    } else {
        let mut l_p = BigInt::from(0);
        let  mi = (nn.clone() - BigInt::from(1)) / BigInt::from(2);
        let mut temp_mi:u32 = 0;
        match mi.to_u32() {
            Some(i) => {
                temp_mi = i;
            },
            None => {
                // println!("APRtest_step4b is unable to convert BigInt to i64");
            }
        }
        let  temp = (qq.pow(temp_mi) % nn.clone() + BigInt::from(1)) % nn.clone();
        if h.clone() % pp.clone() != l_p &&  temp == BigInt::from(0) {
            l_p = BigInt::from(1);
        } else {
            l_p = BigInt::from(0);
        }
        return (true, l_p.clone());
    }

    // return (false, BigInt::from(-1).clone());
}

// fn Pow(p: &mut BigInt,e: &mut BigInt, m: &mut BigInt) ->BigInt {
//     let mut pp = p.clone();
//     let mut ee = e.clone();
//     let mut mm = m.clone();
//     let mut res = BigInt::from(1);
//     println!("{:?},{:?},{:?}",pp, ee, mm);
//     while ee > BigInt::from(0) {
//         if ee.clone() % BigInt::from(2) == BigInt::from(1) {
//             res = res.clone() * pp.clone() % mm.clone();
//         }
//         pp = pp.clone() * pp.clone();
//         ee = ee.clone() / BigInt::from(2);
//         // println!("ka");
//     }
//     return res.clone();
// }

fn aprtest_step4c(p: &mut BigInt,k: &mut BigInt, q: &mut BigInt, n: &mut BigInt) -> (bool, BigInt) {
    let mut pp = p.clone();
    let mut kk = k.clone();
    let mut qq = q.clone();
    let mut nn = n.clone();
    //println!("Step 4c. (p^k, q = {0}^{1}, {2})",pp, kk, qq);

    let  j2q = calc_j(&mut pp, &mut kk, &mut qq);
    // println!("{:?}", J2q);
    // let mut tttt = J2q.clone() * J2q.clone()* qq.clone();
    // println!("{:?}, {:?}", tttt, NN);
    let mut s1 = (j2q.clone() * j2q.clone() * qq.clone()).mmod(&mut nn);
    // println!("{:?}", s1);
    // println!("{:?},{:?}", NN.clone() / BigInt::from(4), NN);
    let  s2 = s1.modpow(&mut (nn.clone() / BigInt::from(4)), &mut nn);
    // println!("{:?}", s2);
    let mut s = JacobiSum::new(&mut pp,&mut kk, &mut qq);
    if nn.clone() % BigInt::from(4) == BigInt::from(1) {
        s = s2.clone();
    } else if nn.clone() % BigInt::from(4) == BigInt::from(3) {
        s = (s2.clone() * j2q.clone() * j2q.clone()).mmod(&mut nn);
    } else {
        println!("Error");
    }
    s = s.mmod(&mut nn);
    // println!("{:?},{:?}", S, NN);
    
    let (exist, h) = s.is_root_of_unity(&mut nn);
    // println!("{:?},{:?}",exist,h);
    if !exist {
        return (false, BigInt::from(0).clone()); // BigInt::from(-1).clone()表示None 
    } else {
        let mut l_p = BigInt::from(0);
        let mut mi = (nn.clone() - BigInt::from(1)) / BigInt::from(2);
        
        // println!("{:?}, {:?}, {:?}", qq, mi, NN);
        let  temp = (qq.modpow(&mut mi, &mut nn) % nn.clone() + BigInt::from(1)) % nn.clone();
        // println!("{:?}", temp);
        if h.clone() % pp.clone() != l_p &&  temp == BigInt::from(0) {
            l_p = BigInt::from(1);
        } else {
            l_p = BigInt::from(0);
        }
        return (true, l_p.clone());
    }

    // return (false, BigInt::from(-1).clone());//BigInt::from(-1).clone()代表None
}

fn aprtest_step4d(q: &mut BigInt, n: &mut BigInt) -> (bool, BigInt) {
   //aprtest_step4d(&mut pp, &mut kk, &mut qq, &mut nn);
    //let  pp = p.clone();
    //let  kk = k.clone();
    //let  qq = q.clone();
    let  nn = n.clone();
    //println!("Step 4d. (p^k, q = {0}^{1}, {2})",pp, kk, qq);
    let  mi = (nn.clone() - BigInt::from(1)) / BigInt::from(2);
    let mut temp_mi:u32 = 0;
    match mi.to_u32() {
        Some(i) => {
            temp_mi = i;
        },
        None => {
            // println!("APRtest_step4b is unable to convert BigInt to i64");
        }
    }
    let fuq = -q.clone();
    let s2q = fuq.pow(temp_mi) % nn.clone();
    if (s2q.clone() -BigInt::from(1)) % nn.clone() != BigInt::from(0) && (s2q.clone() + BigInt::from(1)) % nn.clone() != BigInt::from(0) {
        return (false, BigInt::from(0).clone()); //BigInt::from(-1).clone() represent None
    } else {

        let mut l_p = BigInt::from(0);
        if (s2q.clone() + BigInt::from(1)) % nn.clone() == l_p &&  (nn.clone() - BigInt::from(1)) % BigInt::from(4) == BigInt::from(0) {
            l_p = BigInt::from(1);
        } else {
            l_p = BigInt::from(0);
        }
        return (true, l_p.clone());
    }

}

fn aprtest_step4(p: &mut BigInt,k: &mut BigInt, q: &mut BigInt, n: &mut BigInt) -> (bool, BigInt) {
    let mut pp = p.clone();
    let mut kk = k.clone();
    let mut qq = q.clone();
    let mut nn = n.clone();
    let mut result = false;
    let mut l_p = BigInt::from(0);
    if pp >= BigInt::from(3) {
        (result, l_p) = aprtest_step4a(&mut pp, &mut kk, &mut qq, &mut nn);
        // println!("{:?}, {:?}", result, l_p);
    } else if pp == BigInt::from(2) && kk >= BigInt::from(3) {
        (result, l_p) = aprtest_step4b(&mut pp, &mut kk, &mut qq, &mut nn);
        // println!("{:?}, {:?}", result, l_p);
    } else if pp == BigInt::from(2) && kk == BigInt::from(2) {
        (result, l_p) = aprtest_step4c(&mut pp, &mut kk, &mut qq, &mut nn);
        // println!("{:?}, {:?}", result, l_p);
    } else if pp == BigInt::from(2) && kk == BigInt::from(1) {
        //(result, l_p) = aprtest_step4d(&mut pp, &mut kk, &mut qq, &mut nn);
        (result, l_p) = aprtest_step4d(&mut qq, &mut nn);
    } else {
        println!("error");
    }
    
    if !result {
        // println!("Composite");
    }
    return (result, l_p.clone());
}

fn apr_test(n: &mut BigInt) -> bool {
    let t_list = vec![
        BigInt::from(2),
        BigInt::from(12),
        BigInt::from(60),
        BigInt::from(180),
        BigInt::from(840),
        BigInt::from(1260),
        BigInt::from(1680),
        BigInt::from(2520),
        BigInt::from(5040),
        BigInt::from(15120),
        BigInt::from(55440),
        BigInt::from(110880),
        BigInt::from(720720),
        BigInt::from(1441440),
        BigInt::from(4324320),
        BigInt::from(24504480),
        BigInt::from(73513440)
    ];
    //println!("N={}", N);
    if n <= &mut BigInt::from(3) {
        println!("input should be greater than 3");
        return false
    }
    let mut t_f:i64 = 0; 
    let mut t = BigInt::from(0);
    let mut et = BigInt::from(0);
    let mut q_list = Vec::new();

    for tt in t_list {
        // println!("{}", tt);
        let mut temp = tt.clone();
        let (ett, q_listt) = e(&mut temp);
        et = ett.clone();
        q_list = q_listt.clone();
        // println!("et:{}", ett);
        if n < &mut (ett.clone() * ett.clone()) {
            // println!("{}", ett);
            t_f = 1;
            t = tt;
            break;
        }
    }
    if t_f == 0 {
        println!("t not found");
        return false;
    }

    //println!("t={:?}", t);
    //println!("e(t)={:?} \n {:?}", et, q_list);
    //println!("=== Step 1 ===");
    let t1 = t.clone() * et.clone();
    let g = n.gcd(&t1);
    if g > BigInt::from(1) {
        println!("Composite");
        return false;
    }
    //println!("=== Step 2 ===");
    let mut l: HashMap<BigInt, BigInt> = HashMap::new();
    let mut tt = t.clone();
    let fac_t = prime_factorize(&mut tt);
    //println!("fac_t={:?}", fac_t);
    for (p, _k) in fac_t {
        let mut temp_p:u32 = 0;
        match p.to_u32() {
            Some(i) => {
                temp_p = i;
            },
            None => {
                // println!(" unable to convert BigInt to i64");
            }
        }
        // println!("{:?}", temp_p);
        let t2 = n.pow(temp_p - 1) % (p.clone() * p.clone());
        // println!("{:?}", t2);
        if p >= BigInt::from(3) && t2 != BigInt::from(1) {
            l.insert(p.clone(), BigInt::from(1));
        } else {
            l.insert(p.clone(), BigInt::from(0));
        }
    }
    // println!("l_p={:?}", l);

    let q_listtt = q_list.clone();


    //println!("=== Step 3&4 ===");
    for q in q_list {
        if q == BigInt::from(2) {
            continue;
        }
        let fac =  prime_factorize(&mut (q.clone() - BigInt::from(1)));
        // println!("{:?}, {:?}", q,fac);
        for (p, k) in fac {
            // step4
            let mut pp = p.clone();
            let mut kk = k.clone();
            let mut qq = q.clone();
            let mut nn = n.clone();
            // println!("{:?}, {:?}, {:?}, {:?}",pp, kk, qq,NN);
            let (result, l_p) = aprtest_step4(&mut pp, &mut kk, &mut qq, &mut nn);
            // println!("{:?}, {:?}, {:?}",q, result, l_p);
            if !result {
                println!("Composite");
                return false;
            } else if l_p == BigInt::from(1) {
                // println!("{:?},{:?}",p.clone(), BigInt::from(1));
                l.insert(p.clone(), BigInt::from(1));
            }
        }
    }
    // l.insert(BigInt::from(2), BigInt::from(1));


    //println!("=== Step 5 ===");
    //println!("l_p={:?}", l);    
    for (key, value) in l.iter() {
        if value.clone() == BigInt::from(0) {
            //println!("Try other (p,q). p={:?}", key);
            let mut count:u32 = 0;
            let mut i:u32 = 1;
            let mut found:bool = false;
            while count < 30 {
                let mut q = key.clone() * BigInt::from(i) + BigInt::from(1);
                let mut is_contain:bool = false;
                for z in &q_listtt {
                    if z.clone() == q {
                        is_contain = true;
                        break;
                    }
                }
                if n.clone() % q.clone() != BigInt::from(0) && isprime_slow(&mut q) && (!is_contain) {
                    count += 1;
                    let mut nn = n.clone();
                    let mut pp = key.clone();
                    let mut qq = q.clone() - BigInt::from(1);
                    let k = v(&mut pp, &mut qq);
                    // Step 4
                    let mut kk = BigInt::from(k);
                    
                    let (result, l_p) = aprtest_step4(&mut pp, &mut kk, &mut qq, &mut nn);
    //                 // println!("{:?}, {:?}", result, l_p);
                    if !result {
                        println!("Composite");
                        return false;
                    } else if l_p == BigInt::from(1) {
                        found = true;
                        break;
                    }
                }
                i += 1;
            }
            if !found {
                //println!("error in Step 5");
                return false;
            }
        }
    }

    //println!("=== Step 6 ===");
    let mut r = BigInt::from(1);
    let mut temp_t:i32 = 0;
    match t.to_i32() {
        Some(i) => {
            temp_t = i;
        },
        None => {
            // println!("unable to convert BigInt to i64");
        }
    }

    for _ in 0..temp_t - 1 {
        r = (r.clone() * n.clone()) % et.clone();
        if r.clone() != BigInt::from(1) && r.clone() != n.clone() && r.clone() % n.clone() == BigInt::from(0) {
            //println!("Composite {:?}", r);
            return false;
        }
    }

    println!("Prime");
    return true
}


fn main() {
    // let s = "9";
    // let mut n = BigInt::from_str(s).unwrap();
    // let start = Instant::now();
    // apr_test(&mut n);
    // let duration = start.elapsed();
    // println!("Time: {:?}", duration);

    // Get the input file name from command-line arguments.
    let args: Vec<String> = env::args().collect();
    let input_filename = &args[1];

    // Read the input from the specified file.
    let input_content = fs::read_to_string(input_filename).expect("Failed to read file");

    // Parse the input content into a Vec<BigInt>.
    let input: Vec<BigInt> = input_content
        .trim()
        .split_whitespace()
        .map(|num| BigInt::from_str(num).expect("Failed to parse BigInt"))
        .collect();

    // Find the prime numbers in the input vector.
    let start_time = Instant::now();
    let mut pseudo_primes = 0;
    for n in &input {
        let mut num = n.clone();
        println!("{}", num);
        // If apr_test != true, then n is not prime but pseudo prime.
        if !(apr_test(&mut num) == true){
            pseudo_primes += 1;
        }
    }
    let elapsed_time = start_time.elapsed();

    // Output the count of pseudo primes and the time taken.
    println!("The count of pseudo primes from the previous test: {}", pseudo_primes);
    println!("Time: {:?}", elapsed_time);

}