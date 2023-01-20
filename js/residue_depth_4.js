var select = source.selected.indices;
const wd = source.data;
const vd = other.data;
if (select.length == 1) {
    var resi1 = wd['resi1'][select[0]]
    var resi2 = wd['resi2'][select[0]]
    var counter = 0;
    var idx = 0;
    var tracker = 0;
    for (const i of vd['resi1']) {
        if (i == resi1) {
            counter++;
            idx = tracker;
        }
        tracker++;
    }
    if (counter != 1) {
        counter = 0;
        idx = 0;
        tracker = 0;
        for (const i of vd['resi2']) {
            if (i == resi2) {
                counter++;
                idx = tracker;
            }
        tracker++;
        }
    }
    if (counter == 1) {
        other.selected.indices = [idx];
    }
}