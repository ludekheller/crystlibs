#!/usr/bin/env python3
"""
Script to distribute functions from crystallography_functions_PART*_ACTUAL.py
into the appropriate library files (crystlib, orilib, plotlib, projlib).
"""

import re
import os

# Define function mappings based on their purpose and names
FUNCTION_MAPPINGS = {
    'crystlib': [
        # Lattice vector operations
        'cubic_lattice_vec',
        'monoclinic_lattice_vec',
        'tetragonal_lattice_vec',
        'lattice_vec',
        'reciprocal_basis',
        # Miller index operations
        'uvtw2uvw',
        'uvw2uvtw',
        'hkil2hkl',
        'hkl2hkil',
        # Coordinate transformations
        'xyz2fractional',
        'miller2fractional',
        'xyz2fractional02',
        'normArrayColumns',
        # String conversions
        'vec2string',
        'plane2string',
        'dir2string',
        # Utility functions
        'find_gcd',
        'perpendicular_vector',
        # Tensor operations
        'permut_tensor3',
        'np_permut_tensor3',
        'kronecker',
        'np_kronecker',
        # Symmetry and equivalence
        'symmetry_elements',
        'equivalent_elements',
        # Lattice correspondence
        'B19p_B2_lattice_correspondence',
        'lattice_correspondence',
        'B19p_B2_lattice_correspondence_ini',
        'cubic2tetragonal_lattice_correspondence',
        'Rp_B2_lattice_correspondence',
        'print_correspondence',
        'write_lattice_correspondence',
        # Hexagonal slip systems
        'gensystemsHexIni',
        'gensystemsHex',
        'genallHexSys',
        # Atomic positions and lattice generation
        'generate_lattite_atom_positions',
        'generate_lattice_vectors',
        'generate_lattice_points',
        'generate_lattice_faces',
        'generate_product_lattice_points',
        'generate_product_lattice_faces',
        'generate_plane_vertices',
        # Plane selection
        'select_atomic_plane',
        'select_plane',
        'select_atomic_region',
        'select_crystal_planes',
        # Twinning operations
        'an_between_vecs',
        'habitplane_equation_solution',
        'twinnedhabitplane',
        'twin_equation_solution_ini',
        'twin_equation_solution',
        'get_twinning_plane_points',
        'get_interface2d',
        # Deformation gradients
        'def_gradient_stressfree',
        'def_gradient_stressfree_ini',
        'def_gradient',
        'def_gradient_ini',
        'def_gradient_ini2',
        # NiTi specific
        'niti_twinning',
        'get_twinningdata',
        'get_twinning_dislocation',
        'gen_twinned_lattice_points',
        # File I/O
        'write_txt',
        'read_txt',
        # Utility
        'plane_line_intersection',
        'flipvector',
        'flipvector2negative',
        'vector2miller_ini',
        'vectors2miller',
    ],
    
    'orilib': [
        # Euler angles and rotation matrices
        'np_euler_matrix',
        'np_inverse_euler_matrix',
        # Rodrigues-Frank vectors (axis-angle)
        'ol_g_rtheta_rad',
        'np_ol_g_rtheta_rad',
        'ol_rtheta_g_rad',
        'np_ol_rtheta_g_rad',
        'rotation_from_axis_angle',
        # Rodrigues-Frank vector conversions
        'ol_g_R',
        'np_ol_g_R',
        'ol_R_g',
        'np_ol_R_g',
        'ol_g_R2',
        'np_ol_g_R2',
        'ol_R_g2',
        'np_ol_R_g2',
        # Quaternion operations
        'np_ol_R_q2',
        'np_ol_g_q2',
        'np_ol_q_g',
        # Rotation operations
        'active_rotation',
        'passive_rotation',
        # Euler angle reduction
        'euler_angles_reduction',
    ],
    
    'plotlib': [
        # Mohr circles plotting
        'mohr_circles',
        'plot_mohr_circles',
        'plot_planes_on_mohr_circle',
        'write_mohr_planes',
        'strains_along_13mohrcirle',
        'zero_normal_strains',
        # Lattice plotting 3D
        'plot_lattice3D',
        'plot_latticefaces3D',
        'plot_latticesfaces3D',
        'set_aspect_equal_3d',
        # Lattice plotting 2D
        'plot_lattice2D',
        'plot_lattice_2Dprojection',
        'plot_lattice',
        'plot_lattice_proj',
        'plot_points_proj',
        'plot_lattice_plane',
        'plot_lattice_boundaries',
        # Atomic plane plotting
        'plot_atomic_plane2D',
        'plot_atomic_plane3D',
        'plot_atomlattice2D',
        'plot_cut2D',
        # Stereographic projection plotting
        'plot_planes_on_stereotriangle',
        'plot_planes_on_wulffnet',
        'plot_princip_dir_on_stereotriangle',
        'plot_princip_dir_on_wulffnet',
    ],
    
    'projlib': [
        # Stereographic projections
        'stereoprojection_directions',
        'equalarea_directions',
    ],
}


def extract_function(content, func_name):
    """
    Extract a complete function definition from content.
    Returns the function as a string with its docstring.
    """
    # Pattern to match function definition
    pattern = rf'^def {re.escape(func_name)}\([^)]*\):[^\n]*\n'
    
    match = re.search(pattern, content, re.MULTILINE)
    if not match:
        return None
    
    start_pos = match.start()
    
    # Find the end of the function (next def at same indentation level or end of file)
    # Look for next function definition at column 0
    next_func_pattern = r'\ndef [a-zA-Z_][a-zA-Z0-9_]*\('
    next_match = re.search(next_func_pattern, content[start_pos + len(match.group(0)):])
    
    if next_match:
        end_pos = start_pos + len(match.group(0)) + next_match.start()
    else:
        end_pos = len(content)
    
    function_text = content[start_pos:end_pos]
    
    # Clean up trailing whitespace but keep function structure
    function_text = function_text.rstrip() + '\n\n'
    
    return function_text


def read_file(filepath):
    """Read file content."""
    with open(filepath, 'r', encoding='utf-8') as f:
        return f.read()


def write_file(filepath, content):
    """Write content to file."""
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)


def get_file_header(library_name):
    """Get appropriate header for each library file."""
    headers = {
        'crystlib': '''#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Crystal Structure and Lattice Operations Library

This module provides comprehensive utilities for crystallographic calculations including:
- Lattice vector generation for all crystal systems
- Miller index operations and conversions
- Coordinate transformations
- Symmetry operations
- Lattice correspondence matrices
- Twinning and deformation analysis
- Atomic position generation

Extended with functions from crystallography_functions.py
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

''',
        'orilib': '''#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Orientation Library Module

This module provides comprehensive utilities for crystallographic orientation analysis including:
- Quaternion operations (conversion, multiplication, misorientation calculations)
- Euler angle and rotation matrix conversions
- Rodrigues-Frank vector representations
- Orientation sampling and gridding (HEALPix, Hopf coordinates)
- Crystal symmetry operations
- Axis-angle representations

Extended with functions from crystallography_functions.py
"""

import numpy as np
import math
from numba import njit

''',
        'plotlib': '''#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Crystallographic Plotting Library Module

This module provides comprehensive plotting utilities for crystallographic data including:
- Stereographic projections (equal-area Schmidt net, equal-angle Wulff net)
- Pole figures and inverse pole figures
- Orientation density plots and histograms
- Color maps for orientation data
- Mohr circles for strain analysis
- Lattice visualization in 2D and 3D
- Atomic plane plotting

Extended with functions from crystallography_functions.py
"""

from projlib import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors

''',
        'projlib': '''#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Stereographic Projection Library Module

This module provides utilities for stereographic projections and pole figures including:
- Equal-area (Schmidt) projections
- Equal-angle (Wulff) projections
- Coordinate transformations
- Pole figure generation

Extended with functions from crystallography_functions.py
"""

import numpy as np
import matplotlib.pyplot as plt

'''
    }
    return headers.get(library_name, '')


def main():
    """Main distribution function."""
    
    base_dir = '/mnt/user-data/outputs'
    
    # Read all PART files
    print("Reading PART files...")
    part1_content = read_file(f'{base_dir}/crystallography_functions_PART1_ACTUAL.py')
    part2_content = read_file(f'{base_dir}/crystallography_functions_PART2_ACTUAL.py')
    part3_content = read_file(f'{base_dir}/crystallography_functions_PART3_ACTUAL.py')
    
    # Combine all content
    all_content = part1_content + '\n' + part2_content + '\n' + part3_content
    
    # Read existing library files
    print("Reading existing library files...")
    existing_libs = {}
    for lib in ['crystlib', 'orilib', 'plotlib', 'projlib']:
        filepath = f'{base_dir}/{lib}_commented.py'
        if os.path.exists(filepath):
            existing_libs[lib] = read_file(filepath)
        else:
            existing_libs[lib] = get_file_header(lib)
    
    # Distribute functions
    print("\nDistributing functions...")
    stats = {lib: {'added': 0, 'skipped': 0, 'missing': 0} for lib in FUNCTION_MAPPINGS.keys()}
    
    for library, func_list in FUNCTION_MAPPINGS.items():
        print(f"\n{'='*60}")
        print(f"Processing {library.upper()}")
        print(f"{'='*60}")
        
        new_content = existing_libs[library]
        
        for func_name in func_list:
            func_text = extract_function(all_content, func_name)
            
            if func_text:
                # Check if function already exists in target file
                if f'def {func_name}(' in existing_libs[library]:
                    print(f"  ⚠ SKIPPED: {func_name} (already exists)")
                    stats[library]['skipped'] += 1
                else:
                    new_content += func_text
                    print(f"  ✓ ADDED: {func_name}")
                    stats[library]['added'] += 1
            else:
                print(f"  ✗ MISSING: {func_name} (not found in PART files)")
                stats[library]['missing'] += 1
        
        # Write updated file
        output_path = f'{base_dir}/{library}_EXTENDED.py'
        write_file(output_path, new_content)
        print(f"\n✓ Created: {output_path}")
    
    # Print summary
    print("\n" + "="*60)
    print("DISTRIBUTION SUMMARY")
    print("="*60)
    for library, counts in stats.items():
        print(f"\n{library.upper()}:")
        print(f"  Added: {counts['added']}")
        print(f"  Skipped (already exists): {counts['skipped']}")
        print(f"  Missing (not found): {counts['missing']}")
        print(f"  Total attempted: {len(FUNCTION_MAPPINGS[library])}")
    
    total_added = sum(s['added'] for s in stats.values())
    total_skipped = sum(s['skipped'] for s in stats.values())
    total_missing = sum(s['missing'] for s in stats.values())
    
    print(f"\n{'='*60}")
    print(f"OVERALL TOTALS:")
    print(f"  Total functions added: {total_added}")
    print(f"  Total skipped: {total_skipped}")
    print(f"  Total missing: {total_missing}")
    print(f"{'='*60}")
    
    print("\n✓ Distribution complete!")
    print("\nNew files created:")
    for lib in FUNCTION_MAPPINGS.keys():
        print(f"  - {lib}_EXTENDED.py")


if __name__ == '__main__':
    main()
