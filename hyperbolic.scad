/* [Geometry] */
// Order *p* of rotation around first triangle vertex (0 for ∞)
p = 3; // [0,2,3,4,5,6,7,8]
// Order *q* of rotation around second triangle vertex
q = 2; // [0,2,3,4,5,6,7,8]
// Order *r* of rotation around third triangle vertex. The sum 1/*p*+1/*q*+1/*r* must be < 1.
r = 7; // [0,2,3,4,5,6,7,8]
// Depth of tiling - warning, exponential complexity
depth = 3; // [0,1,2,3,4]


/* [Style] */
// Radius of object, in millimeters
radius = 40;
// Marked side of triangles
marked_side = 2; // [0,1,2]
width = 0.125;

$fs = .5;
$fa = 2;

/* [Hidden] */ // {{{1
// only works for newer OpenSCAD versions...
// function return(x, s="") = echo(str("returning ", s, "=",x)) x;

// TODO:
// - put words-by-length in reverse order
// - allow curvature
// - allow ideal triangles (FIXME)
// - allow outside perimeter (or circle)
// + redo Klein quartic
// + write a generic module list-of-triangles -> style
// determining value type {{{1
// function is_undef(x) = x == undef;
function is_number(x) = (isa_list(x) ? false: x+1 > x);
// function is_string(x) = str(x) == x;
// we could call this is_list, but recent version of OpenSCAD provide
// a primitive is_list function:
function isa_list(x) = (x[0] == x[0]) && (x[0] != undef) && len(x) != undef;

function default(x,y)=x==undef?y:x;

real_precision = 1e-7;
function is_zero(x, e = real_precision, c = 0) =
  x == undef ? undef :
  is_number(x) ? (abs(x) <= e) :
  (c >= len(x)) || (is_zero(x[c], e) && is_zero(x, e, c+1));

// list functions {{{1
function cons(x, tail) = concat([x], tail);
function prefix(list,n)=[ for(i=[0:min(n,len(list)-1)]) list[i] ];
function tail(list, c=1) = [ for(i=[c:1:len(list)-1]) list[i] ];
function sublist(list, start, count) =
  let(n = is_number(count) && count <= len(list)-start ?
    count : len(list)-start)
  [for(i=[0:1:n-1]) list[start+i] ];

function flatten(list, c=0) =
  (c >= len(list) ? [] : concat(list[c], flatten(list, c=c+1)));
// merges two lists while preserving unicity of elements
// (assuming all elements in list1 are already unique)
function merge(list1, list2, c = 0) = c >= len(list2) ?  list1 :
  merge(find(list2[c], list1) ? list1 : concat(list1, [list2[c]]), list2, c+1);

function belongs(z, list, c = 0) =
  (c < len(list)) && (is_zero(list[c]-z) || belongs(z, list, c+1));
function find(x, list, c=0) = c > len(list) ? false :
  list[c] == x ? true : find(x, list, c+1);
function assoc(list, key, c=0) = c >= len(list) ? undef :
  list[c][0] == key ? list[c][1] : assoc(list, key, c+1);

function sum(list, c = 0) = c >= len(list) ? 0 : list[c] + sum(list, c + 1);
function norm2(list, c = 0) = list[c]*list[c] +
  (c < len(list) - 1 ? norm2(list, c + 1) : 0);

function repeat(what, count) = [ for(i=[1:1:count]) what ];
function repeat_list(what, count) = flatten([ for(i=[1:1:count]) what ]);
// function rotate_right(list, k=1) =
//   // rotate_right[a0, a1,...] = [a(n-1), a0, ... ]
//   // we allow shifts to be any integer, even <0 or >= n:
//   let(n = len(list), k=mod(k,n))
//   concat([for(i=[n-k:1:n-1]) list[i]], [for(i=[0:1:n-k-1]) list[i]]);

// Primitive matching engine //{{{1
function is_suffix(list,pattern,c=0,d) =
  let(d = default(d, len(list) - len(pattern)))
  c >= len(pattern) ||
  (list[c+d] == pattern[c] && is_suffix(list, pattern, c+1, d));
function has_suffix_in(list,patterns,c=0) =
  (c < len(patterns)) &&
  (is_suffix(list,patterns[c]) || has_suffix_in(list, patterns, c+1));
function match_tail(list,pattern,x,y) =
/* Checks if the given list tail-matches the pattern;
   if true, returns the position of match start, else return false.

   A pattern is represented as a list, where each element is either
   a scalar (representing itself), or a list (representing the Kleene
   repetition of a pattern). Disjunction is not implemented (but see
   below).

   Parameters x and y are the position in the list and the pattern.
*/
//   echo(str("match_tail:", list, pattern, x, y))
  let(x = default(x, len(list) - 1),
      y = default(y, len(pattern) - 1))
//   echo(str("  setting x=",x,"  and y=",y))
  y < 0 ? x+1 :
  x < 0 ? false :
  isa_list(pattern[y]) ?
//   echo("  is list")
  // x = p q* iff x = p or x = x'q and x' = pq*
  // we take the not-greedy approach:
  let(r = match_tail(list, pattern, x, y-1))
//   echo(str("   returned ", r))
  (r != false ? r :
  let(z = match_tail(list, pattern[y], x))
//   echo(str("   z = ", z))
  (z == false ? z : match_tail(list, pattern, z-1, y))) :
  list[x] == pattern[y] ?
//   echo("  char match, looking at previous position")
  match_tail(list, pattern, x-1, y-1) : false;
  ;
function match_tail_in(list,patterns,c=0) =
/* Checks if the list matches one of the provided patterns.
   If it matches, return a value evaluating as boolean true; else return
   false.
*/
  (c < len(patterns)) &&
  (match_tail(list,patterns[c]) != false || match_tail_in(list, patterns, c+1));

// complex numbers {{{1
/* Complex numbers are represented as a pair [real, imag].
*/
C_i = [0,1];
function is_complex(x) = len(x)==2 && is_number(x[0]) && is_number(x[1]);

// constructor and accessors:
function C_complex (a,b) =
  is_number(a) ? [a,b] :
  [ for(i=[0:len(a)-1]) C_complex(a[i], b[i]) ];
// the accessors C_real and C_imag also work for real numbers
function C_real(z) = is_number(z) ? z :
  (is_complex(z) ? z[0] : [ for (x = z) C_real(x) ]);
function C_imag(z) = is_number(z) ? 0 :
  (is_complex(z) ? z[1] : [ for (x = z) C_imag(x) ]);
function C_conj(z) = is_complex(z) ? [z[0], -z[1]] :[ for (x = z) C_conj(x) ];
// this also works for matrices:
function C_mul(z,w) = C_complex(C_real(z)*C_real(w) - C_imag(z)*C_imag(w),
  C_real(z)*C_imag(w) + C_imag(z)*C_real(w));
function C_muln(z1=1,z2=1,z3=1,z4=1,z5=1,z6=1,z7=1,z8=1) =
  C_mul(C_mul(C_mul(z1,z2),C_mul(z3,z4)),C_mul(C_mul(z5,z6),C_mul(z7,z8)));
// function C_mul(z,w) = [z[0]*w[0]-z[1]*w[1], z[1]*w[0]+z[0]*w[1]];
function C_inv(z) = [z[0], -z[1]]/norm2(z);
function C_div(z,w) = [z[0]*w[0]+z[1]*w[1], z[1]*w[0]-z[0]*w[1]]/norm2(w);
function C_is_real(z) = is_zero(C_imag(z));
function C_powi(z, n) = (n == 0) ? 1 :
  ( (n == 1) ? z :
  let (w=C_powi(z, floor(n/2)), q = C_mul(w, w))
  (n % 2 == 1) ? C_mul(z, q) : q);
function C_homography(M,z) = is_number(M) ? z :
  (is_complex(z)
  ? C_div (C_mul(M[0][0], z) + M[0][1], C_mul(M[1][0], z) + M[1][1])
  : [ for (x=z) C_homography(M, x) ]);
function C_identity(n=2) = [for (i=[0:n-1]) [for (j=[0:n-1])
  [(i==j)?1:0, 0]]];

function C_product(list,start=1,i=0) = i>= len(list) ? start :
  C_mul(list[i], C_product(list, start, i+1));
function C_powers(z,n,start=1,r=[]) =
  // returns (precomputes) the list [z^0, z^1, ..., z^n ]
  len(r) > n ? r : C_powers(z, n, start,
  concat(r, [ r == [] ? start : C_mul(r[len(r)-1], z) ]));

// math functions {{{1
function cosh(x) = (exp(x)+exp(-x))/2;
function sinh(x) = (exp(x)-exp(-x))/2;
function tanh(x) = (exp(2*x)-1)/(exp(2*x)+1);
function acosh(x) = log(x+sqrt(x*x-1));
function asinh(x) = log(x+sqrt(x*x+1));
function atanh(x) = .5*log((1+x)/(1-x));
function gcd2(x,y) = (x < y) ? gcd2(y,x) : (y == 0) ? x : gcd2(y, x%y);
function gcd(x1,x2,x3=0,x4=0) = gcd2(gcd2(x1,x2),gcd2(x3,x4));
function mod(x,y) = x >= 0 ? x%y : y-((-x)%y);

function quadratic_root_nearest(a,b,c, y = 0) =//{{{
// solves quadratic equation a x^2 + b x + c = 0;
// returns solution nearest to y
  a < 0 ? quadratic_root_nearest(-a,-b,-c, y) :
  // here we assume **a > 0 **
  let (D=b*b-4*a*c,
  // we do not take sign in case of y == (-b/(2*a)):
      s = (y >= (-b/(2*a)) ? 1 : -1))
      (-b+s*sqrt(D))/(2*a);//}}}
// plane geometry {{{1
// unit vector parallel to I*v
function normal(v) = let (w = [ -v[1], v[0] ]) w/norm(w);
function angle(v) = let(n = norm(v)) acos(v[0]/n) * (v[1] < 0 ? -1 : 1);
function is_on_circle(c,r,x) = is_zero((x-c)*(x-c)-r*r);
function is_on_line(a,b,x) =
  is_zero((x[0]-a[0])*(b[1]-a[1])-(x[1]-a[1])*(b[0]-a[0]));
function is_on_path(p,x) = len(p)==2 ? is_on_line(p[0],p[1],x) :
  is_on_circle(p[2],p[3],x);

function circle_inter_line(c,r,a,b, y) =//{{{
/* Computes the intersection of the circle with center c and radius r
   and the line through points a and b.

   This returns the intersection point closest to y.
   (a, b, c, y are pairs of real numbers).
*/
  let (u = b-a,
  // the solution is z = a + t u, such that
  // (z-c)^2 = r^2, or (a-c+t u)^2 - r^2 = 0, or
  // u^2 t^2 + 2 (a-c)u t + (a-c)^2 -r^2 = 0
    t = quadratic_root_nearest(norm2(u), 2*(a-c)*u, (a-c)*(a-c)-r*r,
      (y-a)*u/(u*u)))
    a + u * t;//}}}
function circle_inter_circle(c1,r1,c2,r2,y) =//{{{
/* Computes the intersection of the circle with center c1 and radius r1
   and the circle with center c2 and radius r2.

   This returns the intersection point closest to y.
   (c1, c2, y are real numbers).
*/
  let(u = c1-c2,
      t = 1/2-(r1*r1-r2*r2)/2/norm2(u),
      a = t*c1+(1-t)*c2,
      b = a + [-u[1], u[0]])
  circle_inter_line(c1,r1,a,b,y);//}}}
function line_inter_line(a1,b1,a2,b2) =//{{{
/* Computes the intersection of the line through points a1 and b1
   and the line through points a2 and b2.

   (a1, a2, b1, b2 are pairs of real numbers).
*/
 // formulas computed with Pari/GP
  let(d = (-a2[1] + b2[1])*a1[0] + (a2[1] - b2[1])*b1[0]
        + (a2[0] - b2[0])*a1[1] + (-a2[0] + b2[0])*b1[1])
  [ ((-a2[0] + b2[0])*b1[1] + (b2[1]*a2[0] - a2[1]*b2[0]))*a1[0]
  + ((a2[0] - b2[0])*a1[1] + (-b2[1]*a2[0] + a2[1]*b2[0]))*b1[0],
  (-a2[1] + b2[1])*b1[1]*a1[0] + ((a2[1] - b2[1])*a1[1]*b1[0]
  + ((b2[1]*a2[0] - a2[1]*b2[0])*a1[1]
  + (-b2[1]*a2[0] + a2[1]*b2[0])*b1[1]))] / d;//}}}
function homography(M, z) =
  // TODO add a special marker for infinity ([] looks good)
  (M[0][0]*z+M[0][1]) / (M[1][0]*z + M[1][1]);


// Simple substitution in words {{{1
// function word_subs(word, before, after, pos=0, match=0) =
//   // we look for a match starting at position [pos], knowing
//   // that already [match] letters match.
//   (match == 0 && pos + len(before) > len(word)) ? word :
//   match >= len(before) ?
//     concat(sublist(word, 0, pos), after,
//            sublist(word, pos+len(before))) :
//   word[pos+match] == before[match] ?
//     word_subs(word, before, after, pos, match+1) :
//     word_subs(word, before, after, pos+1, 0);
// function word_subs_once(word, rules, start=0) =
//   start >= len(rules) ? word :
//   word_subs_once(word_subs(word, rules[start][0], rules[start][1]),
//     rules, start+1);
// function word_subs_stable(word, rules) =
//   let(new_word = word_subs_once(word, rules))
//   word == new_word ? word : word_subs_stable(new_word, rules);
// Generating words according to some substitution rules {{{2
// function rules_new_words(word_list, alphabet, rules) =
//   flatten([ for(w = word_list) [ for(a = alphabet)
//     word_subs_stable(concat(w, [a]), rules) ] ]);
// function rules_gen_words(alphabet, rules, length=2, previous = [ [] ]) =
//   length == 0 ? previous :
//   merge(previous,
//     rules_gen_words(alphabet, rules, length-1,
//     rules_new_words(previous, alphabet, rules)));
// hyperbolic geometry (Poincaré disk model) {{{1
// Mobius transform such that f(a) = 0
function moebius(a) = [ [[1,0], -a], [-C_conj(a), [1,0]]];
function hyperbolic_reflection(a, b) =//{{{
/* Computes the the hyperbolic reflection ρ around the geodesic [a,b].
   (a, b are complex numbers).

   Returns the matrix M such that ρ(z) = homography(M, conj(z)).
*/
  let (z = (b-a) + C_muln(a,b,C_conj(a-b)),
    w=C_mul(a,C_conj(b)) - C_mul(b, C_conj(a)))
  [[z, w], [C_conj(w), C_conj(z)]];//}}}
function hyperbolic_arc(z, w, delta = 0) =//{{{
/* Computes the hyperbolic arc (= geodesic segment) joining z to w,
   offset by distance delta to the left of the path.
   (z, w are complex numbers; delta is a real number).

   Returns the arc in either the form
   [z1, z2]: a straight line segment
   [z1, z2, c, r]: an arc of a circle (always < 180°).
   (c is the center, r is the radius).
*/
  C_is_real(C_mul(w, C_conj(z))) ?
    let (n = delta*normal(w-z)) [z+n, w+n]
    :
    let (c = C_mul( C_i/2/C_imag(C_mul(z, C_conj(w))),
        w*(1+norm2(z)) - z*(1+norm2(w))),
      r = norm(c-z),
      u = z-c, v = w-c,
      d = delta * ((u[0]*v[1]-v[0]*u[1] >= 0) ? -1 : 1))
    [z*(1+d/r)-c*d/r, w*(1+d/r)-c*d/r, c, r+d];//}}}
function arc_inter(p1,p2, y) =//{{{
/* Computes the intersection of two hyperbolic arcs in the previous form
   (i.e. each arc p1 and p2 is either [z1, z2] or [z1, z2, c, r]).
   Returns the intersection point closest to y (complex number).
*/
  len(p1) == 2
    ? (len(p2) == 2
     ? line_inter_line(p1[0], p1[1], p2[0], p2[1])
     : circle_inter_line(p2[2],p2[3], p1[0],p1[1], y))
    : (len(p2) == 2
     ? circle_inter_line(p1[2], p1[3], p2[0], p2[1], y)
     : circle_inter_circle(p1[2], p1[3], p2[2], p2[3], y));//}}}
function arc_path(p, radius=100, closed=false) =//{{{
/* OBSOLETE: traces the arc.
*/
// p is either [a,b] or [z,w,c,r]
// in the first case, we return just the line segment;
// in the second one, we need to trace a circle segment (the smaller one)
  radius * (len(p) == 2 ? (closed ? p : [p[0]]) :
  // compute determinant of [z-c, w-c] to obtain the sign
  let (z=p[0], w=p[1], c=p[2], r=p[3],
    u = z-c, v = w-c,
    // use $fa (angle), $fs (size) to determine the number of points
    step = max($fa/4, $fs*radius/r),
    // coordinates of vū
    vu = [u[0]*v[0]+u[1]*v[1], u[0]*v[1]-u[1]*v[0]],
    start = angle(u),
    diff = angle(vu),
    n = ceil(abs(diff)/step))
      // we omit the end point (i=n), since we will output closed loops
      // this will avoid repetitions
      [ for (i = [0:n-(closed?0:1)]) let (t = start+(i/n)*diff)
        c + [cos(t), sin(t)]*r ]);//}}}
function standard_incenter(x, a, c) =//{{{
  // returns the incenter of the triangle ABC,
  // where A = [x,0], angle a, C = [0,0], angle c
  let( u=cos(a/2), v=sin(a/2),
    t=(1-x*x)/(2*x*v),
    z = [x,0] + t*[v, u])
  circle_inter_line (z,t, [0,0], [cos(c/2),sin(c/2)], [x,0]);//}}}
function arc_cut(p,a,b) = 
  len(p) == 2 ? [a, b] : [a, b, p[2], p[3]];
function arc_cut_path(p,a,b,radius) = arc_path(arc_cut(p,a,b),radius);
// // breadth-first tree traversal {{{1
// function homography_only_new(list, visited, base_point, i=0) =//{{{
//   // list is a list of [ new matrix, other_info... ]
//   // (currently other_info = word to new matrix + image of center)
//   // returns [new candidates in this list, all visited points]
//   (i >= len(list)) ? [ [], visited ] :
//   let (rest = homography_only_new (list, visited, base_point, i+1),
//     candidate = list[i],
//     test = candidate[2]
// //     test = C_homography(candidate[0], base_point)
//     )
//     belongs(test, rest[1]) ? rest :
//     [ cons(candidate, rest[0]), cons(test, rest[1]) ];//}}}
// function homography_tree (generators, limit, base_point = [.1, .1]) =//{{{
//   // returns [ nodes, visited ], where
//   //  - nodes is a list of [ matrix, generating-path ]
//   //  - visited is a list of already visited points (matching nodes)
//   let(init = [C_identity(2), [], base_point])
//   homography_tree_level (generators, 1, [init], [base_point],
//     [init], base_point, limit);
// function homography_tree_level (generators, depth, nodes,
//   visited, leaves, base_point, limit) =
//   // nodes = all visited nodes, as [ matrix , [ i1, .., in ] ]
//   //   where matrix = step[i1] ... step[in] (product in this order)
//   // visited = all visited points
//   // leaves = active nodes only
//   depth > limit ? nodes :
//   let (next = flatten([ for (x = leaves) [ for (j = [0:len(generators)-1])
//     let (matrix = C_mul(generators[j], x[0]))
//     [ matrix, cons(j, x[1]), C_homography(matrix, base_point) ] ] ]))
//   let (new = homography_only_new (next, visited, base_point))
// //   echo(str("  new[0] = ", [for(n=new[0]) mtr(n)]))
// //   echo(str("  new[1] = ", [for(n=new[1]) n]))
//   // new = [ new nodes, all visited points ]
//   homography_tree_level (generators, depth+1,
//     concat(nodes, new[0]),
//     concat(visited, new[1]), new[0], base_point, limit);//}}}
// function all_edges(list) =//{{{
//   // takes a list of triangles,
//   // and returns all edges of the list
//   // triangle [ path, side, orientation, word ]
//   // edge = [ A, B, orientation, side ]
//   flatten([ for (x = list)
//     let (path = x[0], orientation=x[1])
//     [ for(e = [0:2]) [ path[(e+1)%3], path[(e+2)%3], path[e],
//       e, orientation ] ] ]);
// function edge_find(i, list, j=0) =
//   j == i ? edge_find(i, list, j+1) :
//   j >= len(list) ? false :
// //   echo(str("is_zero: ", list[i][0], list[i][1], list[j][0], list[j][1]))
//   is_zero(list[i][0]-list[j][0]) && is_zero(list[i][1]-list[j][1]) ? j
//   : edge_find(i, list, j+1);
// function edge_match(edges) =
// //   echo(str("edge match on ", len(edges), " edges"))
//   [ for (i=[0:len(edges)-1])
//     concat(edges[i], [edge_find (i, edges)]) ];//}}}
// Triangle groups {{{1
// this is (e^(acosh(x))-1)/(e^(acosh(x)+1))
function hyp_point_at(x) = sqrt((x-1)/(x+1));
function standard_triangle(p,q,r) =//{{{
/* Computes a standard triangle with angles π/p, π/q, π/r.
   Returns a list [A,B,C,P], where A,B,C are the vertices
   and P is an inner point of the triangle (incenter if possible).

   C is the origin point [0,0] if possible.
*/
  p > 0 ? let (a=180/p, b=180/q, c=180/r,
    A = (cos(b)*cos(c)+cos(a))/(sin(b)*sin(c)),
    B = (cos(c)*cos(a)+cos(b))/(sin(c)*sin(a)),
    y = hyp_point_at(B),
    incenter = standard_incenter (y, 180/p, 180/r))
    [ y*[1,0], hyp_point_at(A)*[cos(c),sin(c)], [0,0],
      incenter ]
  : q > 0 ? let (b=180/q, c=180/r,
    A = (cos(b)*cos(c)+1)/(sin(b)*sin(c)),
    x = hyp_point_at (A),
    incenter = standard_incenter (x, 180/q, c))
    [ [1,0], x* [cos(c), sin(c)], [0,0], [incenter[0], -incenter[1]]]
  : r > 0 ? let (c=180/r, t=tan(c/2),
    incenter =(sqrt(t*t-1)-t)*[cos(c/2),sin(c/2)])
    [ [1,0], [cos(c), sin(c)], [0,0], incenter]
  : let (y=sqrt(3)/2)
    [ [1,0], [-.5,y], [-.5,-y], [0,0] ];
  //}}}
// Constructor and accessors//{{{
// a triangle group G is encoded as:
// G[0]: the angles (3 integers)
// G[1]: the basic triangle (3 vertices + 1 incenter)
// G[2]: the rotations (6 matrices: P, Q, R, P^-1, Q^-1, R^-1)
// G[3]: the reflections
function triangle_group(p,q,r) =
  let(T = standard_triangle(p,q,r))
  let(S = [ for (i=[0:2]) hyperbolic_reflection(T[(i+1)%3], T[(i+2)%3]) ])
  let(R = concat(
    [ for (i=[0:2]) C_mul(S[(i+1)%3], C_conj(S[(i+2)%3])) ],
    [ for (i=[0:2]) C_mul(S[(i+2)%3], C_conj(S[(i+1)%3])) ]))
  [[p,q,r], T, R, S];
// infinite_period = 4;
function tgroup_angle(G, i) = let(t = G[0][i%3]) t>0?t:infinite_period;
function tgroup_base_triangle(G) = G[1];
function tgroup_rotation(G) = G[2];
function tgroup_reflection(G) = G[3];
//}}}
function tgroup_neighbours(G, side = undef) =//{{{
  // returns a list [ [ M, e ]*, d ], such that
  // - d is an integer and (M ↦ e) is a group morphism PGL₂(ℂ) → (ℤ/d)
  // if side == undef, we return the 6 immediate neighbours (by rotations
  //    R_i, R_i^(-1));
  // else we return the neighbours for the full polygon by R[side].
  let (R = tgroup_rotation(G))
  (side == undef) ? [ [ for (i=[0:5]) [R[i], 0] ], 1] :
  let (p = (side+1)%3, q=(side+2)%3,
    d=gcd(tgroup_angle(G, p), tgroup_angle(G, q)))
  [ flatten ([ for (k=[0:1:tgroup_angle(G, side)-1])
    let (Rk = C_powi(R[side], k))
    concat(
    // here we use the relation PQR = 1
      [ for (i = [1:1:tgroup_angle(G, p)-2])
        [C_mul (Rk, C_powi (R[p], i)), (i % d) ] ],
      [ for (j = [1:1:tgroup_angle(G, q)-2])
        [C_mul (Rk, C_powi (R[q], j)), (-j % d) ] ]) ]),
  d ];//}}}
// Substitution rules for the triangle group {{{1
function normalize_pattern_rec(pat,r,k=0,i=0) =//{{{
/* Normalizes the pattern R^*k* *pat* by substituting powers of R
   by letters of the alphabet Σ = {R^k P, R^k Q}.
   *r* is the order of R.
   *k* is the power of R multiplying *pat* on the left.
   *i* is the current position in the pattern list (for pseudo-iteration).

   This returns a list of pairs [ q, s ],
   such that R^k pat is equivalent to the union of patterns q R^s.

   The following rules are applied (L = letter, W = word, s = integer):
   R^k L → (R^k L),
   R^k (L W R^s)* -> R^k | (R^k L) (W R^s L)* W R^s
*/
  i>= len(pat) ? [ [[], k] ] : let(x = pat[i])
  isa_list(x) ? // star case
    let(P = normalize_pattern_rec(x, r, 0, 0))
    flatten([ for (p = P)
    // k (y u)* z = k z | k y (u y)* u z
      let(y = p[0], u = p[1], y1 = tail(y))
      let(uy = concat(mod(y[0]+2*u, 2*r), y1),
          ty = concat(mod(y[0]+2*k, 2*r), y1),
          tz = normalize_pattern_rec(pat, r, k, i+1),
          uz = normalize_pattern_rec(pat, r, u, i+1))
      concat(tz, [ for(a = uz)
        [ concat(ty, [ uy ], a[0]), a[1] ] ])
    ]) :
  x >= 2*r ? normalize_pattern_rec(pat, r, k+x, i+1) :
  let(Q = normalize_pattern_rec(pat, r, 0, i+1))
  [ for(q=Q) [ concat([mod(x+2*k, 2*r)], q[0]), q[1]] ];
function normalize_pattern(pat,r) =
/* Normalizes the pattern *pat* by substituting powers of R.
   Returns a list of patterns whose union is equivalent to *pat*.
*/
  let(v = normalize_pattern_rec(pat,r,0,0))
  [ for(x = v) x[0] ];
  //}}}
function tgroup_patterns_all(G) = ////{{{
/* Given a triangle group G, returns the list of all patterns
   that can be reduced according to the relations in G.

   Patterns are encoded as lists whose elements are either:
    - a number in [0, 2r-1], meaning a letter in the alphabet Σ 
      in the ordering 0 = P, 1 = Q, 2 = RP, 3 = RQ…, 2r-1;
    - a number in [2r, 3r-1], meaning a power of
      encoded as 2r+1 = R, … 3r-1 = R^{-1};
    - or a list, representing the Kleene repetition of another pattern.
   For example, the list [ 0, [2, 1], 2r+1] represents P (RP Q)* R.

   This list is later post-processed by other functions:
    - the powers of R are normalized by normalize_pattern_*;
    - this function actually only returns a set of representatives of the
      patterns modulo left-multiplication by R, which is later expanded
      by the patterns_level function.
*/
  let(p=G[0][0], q=G[0][1], r=G[0][2])
  let(a=floor(p/2), b=floor(q/2), c=floor((p+q)/2))
  let(ra=p-2*a, rb=q-2*b, rc=p+q-2*c)
  [p,q] == [0,0] ? [ [ [0,1], 2], [ [1,2], 2] ] :
  p == 0 ? [ [ [0,1], 2], [ [1,2], 2],
    [ repeat(1, b+1), 2-rb ], [ repeat(2, q-b), rb ], ] :
  let(
    W = (p == 2) ?
        concat(repeat(2*r-2, b-1), 3*r-1, repeat(2*r-2, q-b-1), 3*r-1) :
        concat(repeat(0, a-1),       3*r-1, repeat(1, c-a-1),
               repeat(0, p-(c-b)-1), 3*r-1, repeat(1, q-b-1)),
    X = concat(repeat(3, p-a-1), 2*r+1, repeat(2, q-(c-a)-1),
               repeat(3, c-b-1), 2*r+1, repeat(2, b-1)),
    Y = (q == 2) ?
        concat(repeat(2*r-1, a-1), 3*r-1, repeat(2*r-1, p-a-1), 3*r-1) :
        concat(repeat(1, b-1),       repeat(0, c-b-1), 3*r-1,
               repeat(1, q-(c-a)-1), repeat(0, p-a-1), 3*r-1),
    Z = concat(repeat(2, q-b-1), repeat(3, p-(c-b)-1), 2*r+1,
               repeat(2, c-a-1), repeat(3, a-1), 2*r+1))
  concat([
    [ concat(repeat(3, p-a)), ra ],
    [ concat(repeat(3, p-a-1), 2*r+1, repeat(2, q-c+a)), rc ],
    [ concat(repeat(3, p-a-1), 2*r+1, repeat(2, q-c+a-1),
      repeat(3, c-b)), rb ],
    [ concat(X, [X, 2]), 0 ],
    [ concat(repeat(2, q-b)), rb ],
    [ concat(repeat(2, q-b-1), repeat(3, p-c+b)), rc ],
    [ concat(repeat(2, q-b-1), repeat(3, p-c+b-1), 2*r+1,
      repeat(2, c-a)), ra ],
    [ concat(Z, [Z, 1]), 0 ],
    [ [0, W, 1], 2 ],
    [ concat([0, W], repeat(0, a)), 2-ra ],
    [ concat([0, W], repeat(0, a-1), [3*r-1], repeat(1, c-a)), 2-rc ],
    [ concat([0, W], repeat(0, a-1), [3*r-1], repeat(1, c-a-1),
      repeat(0, p-c+b)), 2-rb ],
    [ [1, Y, 2], 2 ],
    [ concat([1, Y], repeat(1, b)), 2-rb ],
    [ concat([1, Y], repeat(1, b-1), repeat(0, c-b)), 2-rc ],
    [ concat([1, Y], repeat(1, b-1), repeat(0, c-b-1), [3*r-1],
      repeat(1, q-c+a)), 2-ra ],
    ],
    p > 2 ? [] : [
    [ [0, W, 0], 2 ],
    [ concat([0, W], repeat(2*r-2, b)), 2-rb ],
    [ concat([0, W], repeat(2*r-2, b-1), [2*r-1]), 2-rb ],
    ],
    q > 2 ? [] : [
    [ [1, Y, 1], 2 ],
    [ concat([1, Y], repeat(2*r-1, a)), 2-ra ],
    [ concat([1, Y], repeat(2*r-1, a-1), [0]), 2-ra ],
    ]
    );
//}}}
function tgroup_patterns_normalized(G) =//{{{
/* Given a triangle group G, returns a set of normalized patterns
   for the reducible words on Σ according to the relations in G.

   See tgroup_patterns_all for the encoding used for patterns.
   This function returns a list of pairs [pattern, level],
   where the level is an integer indicating the number of letters
   by which the corresponding substitution shortens a word.
*/
  let(r = G[0][2])
  let(L = tgroup_patterns_all(G))
  flatten([ for(x = L) [ for(y = normalize_pattern(x[0], r))
    [ y, x[1] ] ] ]);//}}}
function patterns_level(L, l, r) =//{{{
/* Given a list L of pairs [rep, level]
   where rep is a representative of a class of patterns modulo
   left-multiplication by R (see tgroup_patterns_all),
   and the order r of generator R,
   returns the list of patterns where level >= l
   (obtained by left-multiplication by all powers of R).
*/
  flatten([ for(x = L) x[1] < l ? [] : let(y = x[0])
    [ for(i = [0:1:r-1]) concat(mod(y[0]+2*i, 2*r), tail(y)) ] ]);//}}}
function tgroup_tiles(S, r, length=2) =//{{{
/* Given a list S of avoided patterns and the order r of R,
   returns the set of all words avoiding all patterns in R,
   grouped by word length:
   [ [words of length 0], [words of length 1], …].
*/
  (length == 0) ? [[[]]] : // the only word of length 0
  let(T = tgroup_tiles(S, r, length-1))
  let(W = T[len(T)-1])
  let(new = flatten([ for(g = [0:2*r-1]) flatten([ for(w = W)
    let(w1 = concat(w, [g]))
    match_tail_in(w1, S) ? [] : [w1]
  ]) ]))
  concat(T, [ new ]);//}}}
function tgroup_borders(words, patterns, r) =//{{{
/* Given a set of words for tiles (at a distance n),
   and a set of patterns for borders,
   returns a description of the outer border.

   *r* is the order of R.
   The border is returned as a list of pairs [tile, side],
   where tile is an index in the list of words (provided),
   and side is the index of the generator.
*/
  flatten([ for(i = [0:len(words)-1]) let(w=words[i])
    flatten([for (g = [0:2*r-1])
    let(w1 = concat(w, [g]))
    match_tail_in(w1, patterns) ? [] : [[i, g]]
  ]) ])
;//}}}
// Tessellation data {{{1
function triangle_group_tessellation_data(geometry,l=2)=//{{{
/* Given the parameters [p,q,r] for a triangle group
   and a distance parameter l,
   returns the tessellation data, encoded as a list in the following way:
   [0]: [list of even tiles, list of odd tiles],
     (where a tile is a triple of complex numbers)
   [1]: [list of edges of type 0, list of type 1, list of type 2,
         list of borders]
     (where an edge is a pair of complex numbers)
   [2]: [list of labels]
     (where a label is a pair [string, position = complex number]
*/
  let(G = triangle_group(geometry[0], geometry[1], geometry[2]),
  p=G[0][0], q=G[0][1], r=G[0][2],
  R = tgroup_rotation(G),
  Pi = C_powers(R[0], p-2),
  Qj = C_powers(R[1], q-2),
  Rk = C_powers(R[2], r),
  patterns = tgroup_patterns_normalized(G),
//   a01 = echo([ for(x = patterns) str(x, "\n") ]) 0,
  patterns_tiles = patterns_level (patterns, 0, r),
  patterns_borders = patterns_level (patterns, 1, r),
//   a02 = echo([ for(x = patterns_tiles) str("tile: ", pat_str(x), "\n") ]) 0,
  // compute words:
  words_by_length = tgroup_tiles(patterns_tiles, r, l),
//   words_by_length1 = tgroup_tiles(patterns_tiles, r, l),
//   words_by_length = ([ for(y=(words_by_length1)) flatten([
//     for(x=y) len(x) < 1 || x[0] < 2 ? [x] : []]) ]),
  words = flatten(words_by_length),
  leaves = words_by_length[l],
  first_leaf = sum([ for(i = [0:1:l-1]) len(words_by_length[i]) ]),
  borders = tgroup_borders(leaves, patterns_borders, r),
  // the 2r generators in order:
  generator = flatten([ for(k=[0:1:r-1])
    [ C_mul(Rk[k%r], R[0]), C_mul(Rk[k], R[1]) ] ]),
  tiles = [ for(w = words)
    C_product([ for(g = w) generator[g] ], C_identity(2)) ],
  tilesR = flatten([ for(x = tiles)  [ for(k = [0:1:r-1]) C_mul(x, Rk[k]) ] ]),

  // the paths to draw:
  triangle = tgroup_base_triangle(G),
  odd_triangle = [triangle[1], triangle[2], C_homography(R[2], triangle[0])],
  even_triangles = [ for(y = tilesR)
    [ for(i=[0:2]) C_homography(y, triangle[i]) ] ],
  odd_triangles = [ for(t = even_triangles)
    [ for(v=t) C_homography(tgroup_reflection(G)[0], C_conj(v)) ] ],
  poly = flatten([ for(R=Rk)
    [ for(i=[0:1]) C_homography(R, triangle[i]) ] ]),
  _=0)
//   WE; U=
  [ // [0][0], [0][1] triangles:
    [ even_triangles, odd_triangles, ],
    // [1][0],[1][1], [1][2] edges:
  concat([ for(i = [0:2])
      [ for(t = even_triangles) [t[mod(i+1,3)], t[mod(i+2,3)]] ] ],
    [ [ for(y = borders) let(x = tiles[first_leaf+y[0]], k = y[1])
    // [1][3] borders:
    [ for(i = [0,1]) C_homography(x, poly[mod(k+i-1,2*r)]) ] ] ]),
  // [2] labels:
  [ for(i = [0:len(words)-1])
    [ words[i], C_homography(tiles[i], triangle[2]) ] ],
  ];//}}}
function klein_quartic_data() =//{{{
/* Returns the same data as triangle_group_tessellation_data,
   for the Klein quartic.

   This data is a constant, but it is quite simple to compute once we
   have all the needed functions.
*/
  let(G = triangle_group(2,3,7),
  Q = tgroup_rotation(G)[1],
  R = C_powers(tgroup_rotation(G)[2], 13),
  tiles = concat( [C_identity(2)],
    [ for(i = [0:6]) C_mul(Q,R[i]) ],
    [ for(i = [0:6]) C_muln(Q,R[6],Q,R[i]) ],
    [ for(i = [0:2]) C_muln(Q,R[5],Q,R[i]) ],
    [ for(i = [6:9]) C_muln(Q,R[4],Q,R[i]) ],
    [ C_muln(Q,R[6],Q,R[4],Q,R[1]) ],
    [ C_muln(Q,R[5],Q,R[3],Q,R[1]) ],
  []),
  B = C_homography(C_muln(Q,R[5],Q,R[3],Q,R[1]), [0,0]),
  pi7 = 180/7,
  border = [ for(i = [0:14]) B[0] * [ cos(i*pi7), sin(i*pi7) ] ],
  tilesR = flatten([for(i=[0:1:6]) [ for(x=tiles) C_mul(R[i], x) ] ]),
//     tilesR = flatten([ for(Ri = R) [ for(x = tiles) C_mul(Ri,x) ] ]),
  triangle = tgroup_base_triangle(G),
  even_triangles = [ for(y = tilesR)
    [ for(i=[0:2]) C_homography(y, triangle[i]) ] ],
  odd_triangles = [ for(t = even_triangles)
    [ for(v=t) C_homography(tgroup_reflection(G)[0], C_conj(v)) ] ],
  _=undef)
  [
    [ even_triangles, odd_triangles ],
    concat([ for(i = [0:2])
      [ for(t = even_triangles) [t[mod(i+1,3)], t[mod(i+2,3)]] ] ],
      [[border]]),
    [],
  ];//}}}
// modules for drawing hyperbolic arcs/triangles //{{{1
module hyperbolic_path_draw(x, radius = 100, thickness=0.125) {//{{{
/* Draws the hyperbolic path x (list of complex numbers)
   as a 2-d module, where the Poincaré disk is scaled to given radius.

   This draws len(x)-1 hyperbolic geodesic segments as an open path.
*/
  for(i = [0:len(x)-2]) {
    p = hyperbolic_arc(x[i], x[i+1]);
    if (len(p) == 2)
      hyperbolic_path_draw_line(p[0], p[1], radius, thickness);
    else
      hyperbolic_path_draw_circle (p[0], p[1], p[2], p[3], radius, thickness);
    translate(radius*p[0]) circle(thickness/2);
    translate(radius*p[1]) circle(thickness/2);
  }
}
module hyperbolic_path_draw_line(a, b, radius, thickness) {
  v = thickness/2 * normal(b-a);
  polygon([radius*a+v, radius*b+v, radius*b-v, radius*a-v]);
}
module hyperbolic_path_draw_circle(a, b, c, r, radius, thickness) {
  translate(radius*c) intersection() {
    difference() {
      circle(radius*r+thickness/2);
      circle(radius*r-thickness/2);
    }
    polygon(radius*[[0,0], 2*(a-c), 2*(a+b-2*c), 2*(b-c)]);
  }
}//}}}
module hyperbolic_path_fill(p, orientation, radius=100) {//{{{
/* Fills the hyperbolic path x (list of complex numbers)
   as a 2-dimensional module, where the Poincaré disk is scaled to the
   given radius. The orientation of the path (sign) must be provided.

   This fills the closed path contained by len(x) geodesic segments.
*/
  n = len(p);
  intersection_for(e=[0:n-1]) {
    hyperbolic_left_of(hyperbolic_arc(p[e%n], p[(e+orientation+n)%n]), radius);
  }
}
module hyperbolic_left_of(p, radius, eps=.001) {
  if (len(p) == 2) hyperbolic_left_of_line (p[0], p[1], radius, eps);
  else hyperbolic_left_of_circle(p[0], p[1], p[2], p[3], radius, eps);
}
module hyperbolic_left_of_line(a, b, radius,eps=0) {
  u = (b-a) / norm(b-a) * 3 * radius; // has norm 3*radius
  v = [-u[1], u[0]]; // vector to the left, norm = 3*radius
  off = v*eps/(3*radius);
  intersection() {
    circle(radius);
    polygon([radius*a-u-off, radius*a+u-off, radius*a+u+v-off, radius*a-u+v-off]);
  }
}
module hyperbolic_left_of_circle(a,b,c,r,radius,eps=0) {
  u = c-a; v = c-b;
  if (u[0]*v[1] - u[1]*v[0] > 0) {
    translate(c*radius) circle(r*radius+eps);
  } else {
    difference() {
      circle(radius);
      translate(c*radius) circle(r*radius-eps);
    }
  }
}//}}}
// tessellation modules {{{1
/* These are the basic modules that convert tessellation data into
   three-dimensional objects.

   The first parameter of these modules is the list of things to draw;
   this can be e.g. a list of cells (for hyperbolic_fill)
   or a list of edges (for hyperbolic_draw).
   Appropriate lists are extracted as members of a tessellation data
   structure; see below for examples.

   The next parameters specify the geometry of the modules:
   (e.g. “size” is the radius given to the Poincaré disc).
*/
module hyperbolic_fill(L, orientation=1, size) {//{{{
  color("lightgreen") for(tile = L)
//     linear_extrude(height=height)
    hyperbolic_path_fill(tile, orientation, size);
}//}}}
module hyperbolic_draw(L, size, width=.8) {//{{{
  color("red") for(edge = L)
    // linear_extrude(height=height)
    hyperbolic_path_draw(edge, size, width);
}//}}}
module hyperbolic_text(L, size, txtsize=3) {//{{{
color("black") for(label = L)
//   linear_extrude(height=height)
  translate(size*label[1])
  // pat_str
  text(str((label[0])),size=txtsize);
}//}}}
/*
   These are the higher-level modules converting tessellation data
   into a three-dimensional object.

   The first parameter is a tessellation data structure,
   from which members are extracted and passed to the lower-level modules
   (hyperbolic_draw etc.)

   The next parameters give the geometry of the resulting objects.
*/
module tessellation_simple(T, size = 100) {//{{{
 // T has the form [ tiles, edges, labels ]
 color("lightgreen") linear_extrude(2) hyperbolic_fill(T[0][0], 1, size);
 color("red") linear_extrude(3) hyperbolic_draw(T[1][marked_side], size);
 color("black") linear_extrude(3) hyperbolic_draw(T[1][3], size);
//  hyperbolic_text(T[2], size);
}//}}}
module tessellation_flat(T, size=100) {//{{{
 //hyperbolic_fill(T[0][0], 1, size);
 hyperbolic_draw(T[1][marked_side], size, width);
 //hyperbolic_draw(T[1][3], size, width);
}//}}}
module tessellation_cookie(T, size=50) {//{{{
  linear_extrude(2) hyperbolic_fill(T[0][0], 1, size);
  linear_extrude(3) hyperbolic_draw(T[1][marked_side], size);
  linear_extrude(3) hyperbolic_draw(T[1][3], size);
  // borders
  translate([2*size+10, 0]) {
  color("blue") linear_extrude(1.5)
    hyperbolic_draw(T[1][3], size+1.8, width= 2);
  color("lightblue") linear_extrude(8)
    hyperbolic_draw(T[1][3], size+1.8, width=.8);
  }
}//}}}
module tessellation_bowl(T, size=50) {//{{{
  // Intended to be extruded to a sphere, but see below:
  linear_extrude(2) hyperbolic_fill(T[0][0], 1, size);
  linear_extrude(3) hyperbolic_draw(T[1][marked_side], size);
  linear_extrude(4) hyperbolic_draw(T[1][3], size);
}//}}}
module curved_extrude(r,height=1) {//{{{
/* This is an attempt to do a non-flat module.
   For now, this is beyond the capacities of OpenSCAD.
   In practice, we can still generate a flat STL with OpenSCAD,
   convert this (with admesh) to a text format,
   and then compute orthographic projection on a sphere with another
   program.

   I leave however the code here as an example of how to compute
   in OpenSCAD with modules (where (x,y,z) coordinates are
   inaccessible to the SCAD file).
*/
  if(r == 0) {
    linear_extrude(height) children();
  } else {
    intersection() {
      linear_extrude(height=r) children();
      translate([0,0,r]) difference() {
        sphere(r);
        sphere(r-height);
      }
    }
  }
}//}}}

// {{{1
// }}}1
T = triangle_group_tessellation_data([p,q,r],depth);
// uncomment this line to use the Klein quartic instead:
// T = klein_quartic_data();
//tessellation_cookie(T, 150);
//scale([1.5, 1])
baseRadius = 4.5;
radiusScales = [0.75, 0.5];
xOffsetRatios = [2.58, 3.58];
radiusFudgeFactor =  0.95;
fudgedRadius = baseRadius * radiusFudgeFactor;
offset1 = [-fudgedRadius * xOffsetRatios[0], -fudgedRadius * radiusScales[0]];
offset2 = [-fudgedRadius * xOffsetRatios[1], 0];
/*mirror([1, 0])
intersection(){
    polygon([[-18, 3.5], [18, 3.5], [18,-3.5], [-18, -3.5]]);
    union(){
        //scale ([2, 1])
        translate([-fudgedRadius, 0])
        tessellation_flat(T, baseRadius);
        translate(offset1)
        tessellation_flat(T, baseRadius * radiusScales[0]);
        translate([offset1[0], - offset1[1]])
        tessellation_flat(T, baseRadius * radiusScales[0]);
        translate(offset2)
        tessellation_flat(T, baseRadius * radiusScales[1]);

        //translate([12, 0])
        //tessellation_flat(T, 5);
        //translate([-12, 0])
        //tessellation_flat(T, 5);
    }
};
 intersection(){
    polygon([[-18, 3.5], [18, 3.5], [18,-3.5], [-18, -3.5]]);
    union(){
        //scale ([2, 1])
        translate([-fudgedRadius, 0]) */
        tessellation_flat(T, baseRadius);
        /*translate(offset1)
        tessellation_flat(T, baseRadius * radiusScales[0]);
        translate([offset1[0], - offset1[1]])
        tessellation_flat(T, baseRadius * radiusScales[0]);
        translate(offset2)
        tessellation_flat(T, baseRadius * radiusScales[1]);

        //translate([12, 0])
        //tessellation_flat(T, 5);
        //translate([-12, 0])
        //tessellation_flat(T, 5);
    }
};*/


// vim: fdm=marker:
