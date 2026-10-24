<?php
require __DIR__ . '/../../vendor/autoload.php';
use Rindow\Math\Matrix\MatrixOperator;
use Interop\Polite\Math\Matrix\NDArray;
use function Rindow\Math\Matrix\C;

$mo = new MatrixOperator();
$la = $mo->la();

// 2qubit分のBell状態シミュレーション
$I2 = [[1,0],[0,1]];
function kron($A,$B) {
    $ra=count($A); $ca=count($A[0]); $rb=count($B); $cb=count($B[0]);
    $out=[];
    for($i=0;$i<$ra*$rb;$i++){ $out[$i]=array_fill(0,$ca*$cb,0); }
    for($i=0;$i<$ra;$i++) for($j=0;$j<$ca;$j++) for($k=0;$k<$rb;$k++) for($l=0;$l<$cb;$l++) {
        $out[$i*$rb+$k][$j*$cb+$l] = $A[$i][$j]*$B[$k][$l];
    }
    return $out;
}
$inv = 1/sqrt(2);
$H = [[$inv,$inv],[$inv,-$inv]];
$H_full = kron($H,$I2); // qubit0(MSB)にH、qubit1にI

function toComplexMatrix($m, $la) {
    $rows=[];
    foreach($m as $row){
        $r=[];
        foreach($row as $v){ $r[]=C($v,0); }
        $rows[]=$r;
    }
    return $la->array($rows, dtype: NDArray::complex64);
}
$Hf = toComplexMatrix($H_full,$la);

$psi0 = $la->array([C(1,0),C(0,0),C(0,0),C(0,0)], dtype: NDArray::complex64); // |00>

$psi1 = $la->gemv($Hf, $psi0);
echo "After H(q0): ".$la->toString($psi1)."\n";

// CNOT (control=q0(MSB), target=q1)
$CNOT = [
    [1,0,0,0],
    [0,1,0,0],
    [0,0,0,1],
    [0,0,1,0],
];
$Cf = toComplexMatrix($CNOT,$mo);
$psi2 = $la->gemv($Cf, $psi1);
echo "After CNOT: ".$la->toString($psi2)."\n";

// |amplitude|^2 は la->abs() が complex64 未対応だったため、
// ArrayAccessで real/imag を直接取り出して手計算
//$probs=[];
//for($i=0;$i<4;$i++){
//    $c = $psi2[$i];
//    $probs[] = $c->real*$c->real + $c->imag*$c->imag;
//}
//print_r($probs);

// 確率分布からのショットサンプリング(測定)
//$probsArr = $la->array($probs, dtype: NDArray::float64);
$probsArr = $la->scal(2.0, $la->log($la->abs($psi2)));
echo $la->toString($probsArr)."\n";
$samples = $la->randomCategorical($probsArr, numSamples: 1024);
$counts = $la->bincount($samples,4);
$labels = ['00','01','10','11'];
foreach($counts as $i=>$c){
    echo $labels[$i].": ".$c."\n";
}

// After H(q0): [0.70710678118655+0i,0+0i,0.70710678118655+0i,0+0i]
// After CNOT: [0.70710678118655+0i,0+0i,0+0i,0.70710678118655+0i]
// Array
// (
//     [0] => 0.5
//     [1] => 0
//     [2] => 0
//     [3] => 0.5
// )
// 00: 524
// 01: 0
// 10: 0
// 11: 500
