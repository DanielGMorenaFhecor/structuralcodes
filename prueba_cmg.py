import structuralcodes.codes as c

# import structuralcodes.material.concrete as conc
import structuralcodes.codes.ec2_2004 as ec2

minorma = c._use_design_code('mc2010')

_fck = 25
print('fck = {} MPa'.format(_fck))
print('fcm = {} MPa'.format(minorma.fcm(_fck)))
print('fctm = {} MPa'.format(minorma.fctm(_fck)))
print('Gf = {} MPa'.format(minorma.Gf(_fck)))


# conc.create_concrete(30, 'mihormigon', design_code='fib Model Code 2010')


print(ec2.w_max("XC2", "qp"))
print(ec2._crack_control.w_max("XC2", "qp"))
print(
    'As_min = {} mm2'.format(
        ec2._crack_control.As_min(40000, 500 / 1.15, 2.56, 1, 1)
    )
)
print(
    'As_min = {} mm2'.format(
        ec2._crack_control.As_min(
            40000,
            500 / 1.15,
            minorma.fctm(_fck),
            ec2._crack_control.k_crack_min_steel_area(400),
            ec2._crack_control.kc_crack_min_steel_area_rectangular(
                400, 400, minorma.fctm(_fck), 0
            ),
        )
    )
)
