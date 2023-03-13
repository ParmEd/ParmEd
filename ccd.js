/* 
to run node v16+ required

npm install molstar
node ccd.js components.cif > components.json
*/

const { CIF } = require('molstar/lib/commonjs/mol-io/reader/cif');
const fs = require('fs');

function isValidComponent(block) {
    const { chem_comp } = block.categories;
    const type = chem_comp.getField('type')?.str(0).toLowerCase() ?? '';

    return type.endsWith('linking')
        && (type.includes('dna') || type.includes('rna') || type.includes('peptide'));
}

function formatComponent(block) {
    const { chem_comp, chem_comp_atom, chem_comp_bond } = block.categories;

    // Atoms info
    const atoms = [];
    for (let i = 0; i < chem_comp_atom.rowCount; i++) {
        atoms.push({
            name: chem_comp_atom.getField('atom_id').str(i),
            symbol: chem_comp_atom.getField('type_symbol').str(i),
            formal_charge: chem_comp_atom.getField('charge').int(i),
            leaving: chem_comp_atom.getField('pdbx_leaving_atom_flag').str(i) === 'Y',
            aromatic: chem_comp_atom.getField('pdbx_aromatic_flag').str(i) === 'Y',
        });
        // use .float(i) if float prop is needed
    }

    // Bonds info
    const bonds = [];
    if (chem_comp_bond) {
        for (let i = 0; i < chem_comp_bond.rowCount; i++) {
            let bond = {
                name_a: chem_comp_bond.getField('atom_id_1').str(i),
                name_b: chem_comp_bond.getField('atom_id_2').str(i),
                order: chem_comp_bond.getField('value_order').str(i),
                aromatic: chem_comp_bond.getField('pdbx_aromatic_flag').str(i) === 'Y',
                stereo: chem_comp_bond.getField('pdbx_stereo_config').str(i) === 'Y',
            }
            bonds.push(bond);
        }
    }

    return {
        name: chem_comp.getField('id')?.str(0),
        one_letter_code: chem_comp.getField('one_letter_code')?.str(0),
        full_name: chem_comp.getField('name')?.str(0),
        type: chem_comp.getField('type')?.str(0),
        formal_charge: chem_comp.getField('pdbx_formal_charge')?.int(0),
        atoms,
        bonds,
    };
}

async function run() {
    const data = fs.readFileSync(process.argv[2], 'utf-8');
    const parsed = await CIF.parse(data).run();
    if (parsed.isError) {
        console.error(parsed.message);
        return;
    }

    const components = parsed.result.blocks
        .filter(isValidComponent)
        .map(formatComponent);

    console.log(JSON.stringify(components, null, 2));
}

run();