var select = wild.selected.indices;
const wd = wild.data;
const vd = variant.data;
table.data['wild'] = []
table.data['variant'] = []
for (const j of rows) {
    var total_wd = 0;
    var total_vd = 0;
    for (const i of select) {
        total_wd += wd[j][i];
        total_vd += vd[j][i];
    }
    table.data['wild'].push(total_wd.toFixed(3))
    table.data['variant'].push(total_vd.toFixed(3))
}
table.change.emit();