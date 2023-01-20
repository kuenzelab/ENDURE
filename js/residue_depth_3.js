var j = 0;
for (var i of source.data['label']) {
    if (i == 'Conserved') {
        if (this.active.includes(0)) {
            source.data['color'][j][3] = 1;
        } else {
            source.data['color'][j][3] = 0;
        }
    }
    if (i == 'Better Energy') {
        if (this.active.includes(1)) {
            source.data['color'][j][3] = 1;
        } else {
            source.data['color'][j][3] = 0;
        }
    }
    if (i == 'Worse Energy') {
        if (this.active.includes(2)) {
            source.data['color'][j][3] = 1;
        } else {
            source.data['color'][j][3] = 0;
        }
    }
    j = j + 1;
}
source.change.emit();