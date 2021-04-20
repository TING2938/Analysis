function num(n) {
    var result = new Array();
        result[0] = 1;
        result[1] = 1;
        for  (var i = 2; i < n; ++i) 
            result[i] = result[i-1] + result[i-2]; 
    return result;
}

for (i of num(100)){
    console.log(i);
}
