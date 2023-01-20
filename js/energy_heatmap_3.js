var select = diff.selected.indices;
const wd = diff.data;
table.data['diff'] = []
for (const j of rows) {
    var total_wd = 0;
    for (const i of select) {
        total_wd += wd[j][i];
    }
    table.data['diff'].push(total_wd.toFixed(3))
}
table.change.emit();