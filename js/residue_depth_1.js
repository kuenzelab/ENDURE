var select = source.selected.indices;
const scd = sc.data;
const ocd = oc.data;
if (select.length == 1) {
    var j = select[0] + 1;
    o_s.selected.indices = [j-1];
    scd['resi1'] = [];
    scd['resi2'] = [];
    scd['total'] = [];
    for (const i of s_e[j]) {
        scd['resi1'].push(i[0]);
        scd['resi2'].push(i[1]);
        scd['total'].push(i[2]);
    }
    ocd['resi1'] = [];
    ocd['resi2'] = [];
    ocd['total'] = [];
    for (const i of o_e[j]) {
        ocd['resi1'].push(i[0]);
        ocd['resi2'].push(i[1]);
        ocd['total'].push(i[2]);
    }
    sc.change.emit();
    oc.change.emit();
}
if (select.length == 0) {
    o_s.selected.indices = [];
}