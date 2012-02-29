import sys
import os
sys.path.insert(0,
os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from pyvows import Vows, expect
from pygamess import Gamess
import openbabel as ob
obc = ob.OBConversion()
obc.SetInFormat("mol")
mol = ob.OBMol()
obc.ReadFile(mol, "examples/ethane.mol")
 
@Vows.batch
class PyGamess(Vows.Context):
    class Gamess(Vows.Context):
        def topic(self):
            return Gamess()
 
        def should_be_a_gamess(self, topic):
            expect(topic).to_be_instance_of(Gamess)

        def should_have_a_tempdir(self, topic):
            #expect(topic.tempdir).to_be_str()
            expect(isinstance(topic.tempdir, str)).to_be_true()

        def should_have_a_debug_flag(self, topic):
            expect(topic.debug).to_be_false()

        def should_have_a_path_for_gamess(self, topic):
            expect(topic.gamess).Not.to_be_like("rungms not found")

        def should_have_a_contrl(self, topic):
            expect(isinstance(topic.contrl, dict)).to_be_true()
        def should_have_a_basis(self, topic):
            expect(isinstance(topic.basis, dict)).to_be_true()
        def should_have_a_statpt(self, topic):
            expect(isinstance(topic.statpt, dict)).to_be_true()
        def should_have_a_system(self, topic):
            expect(isinstance(topic.system, dict)).to_be_true()
        def should_have_a_cis(self, topic):
            expect(isinstance(topic.cis, dict)).to_be_true()

        def can_print_section(self, topic):
            expect(topic.print_section).to_be_a_function()
        def print_contrl_should_return_a_text(self, topic):
            expect(topic.print_section('contrl')).to_be_like(' $contrl runtyp=energy scftyp=rhf  $end\n')
        def print_basis_should_return_a_text(self, topic):
            expect(topic.print_section('basis')).to_be_like(' $basis gbasis=sto ngauss=3 $end\n')
        def print_system_should_return_a_text(self, topic):
            expect(topic.print_section('system')).to_be_like(' $system mwords=30  $end\n')
        def print_statpt_should_return_a_text(self, topic):
            expect(topic.print_section('statpt')).to_be_like(' $statpt opttol=0.0001 nstep=20  $end\n')
        def print_cis_should_return_a_text(self, topic):
            expect(topic.print_section('cis')).to_be_like(' $cis nstate=1  $end\n')
        def print_header_should_return_a_text(self, topic):
            expect(topic.print_header()).to_be_like(' $contrl runtyp=energy scftyp=rhf mult=1  $end\n $basis gbasis=sto ngauss=3  $end\n $system mwords=30  $end\n')

        def can_print_gamin(self, topic):
            expect(topic.gamess_input).to_be_a_function()
        def print_gamin_should_return_a_input(self, topic):
            expect(topic.gamess_input(mol)).to_be_like(""" $contrl runtyp=energy scftyp=rhf mult=1  $end
 $basis gbasis=sto ngauss=3 $end
 $system mwords=30  $end
 $DATA
6324
C1
C      6.0     -0.7560000000    0.0000000000    0.0000000000 
C      6.0      0.7560000000    0.0000000000    0.0000000000 
H      1.0     -1.1404000000    0.6586000000    0.7845000000 
H      1.0     -1.1404000000    0.3501000000   -0.9626000000 
H      1.0     -1.1405000000   -1.0087000000    0.1781000000 
H      1.0      1.1404000000   -0.3501000000    0.9626000000 
H      1.0      1.1405000000    1.0087000000   -0.1781000000 
H      1.0      1.1404000000   -0.6586000000   -0.7845000000 
 $END\n\n\n""")
        class WhenChangedABasis(Vows.Context):
            def topic(self, gam):
                return gam
            def should_be_changed_a_basis(self, topic):
                topic.basis_type('am1')
                expect(topic.basis).to_be_like({'gbasis':'am1'})

        class WhenChangedARunType(Vows.Context):
            def topic(self, gam):
                return gam
            def should_be_changed_a_runtyp(self, topic):
                topic.run_type('optimize')
                expect(topic.contrl['runtyp']).to_be_like('optimize')

        class WhenChangedAScf(Vows.Context):
            def topic(self, gam):
                return gam
            def should_be_changed_a_scf(self, topic):
                topic.scf_type('UHF')
                expect(topic.contrl['scftyp']).to_be_like('UHF')

        class AfterRunning(Vows.Context):
            def topic(self, gam):
                return gam.run(mol)
            def should_be_new_mol(self, topic):
                expect(topic.GetEnergy()).to_be_like(-78.30530748)

