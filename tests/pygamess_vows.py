from pyvows import Vows, expect

from pygamess import Gamess
 
@Vows.batch
class PyGamess(Vows.Context):
    class Gamess(Vows.Context):
        def topic(self):
            return Gamess()
 
        def is_a_gamess(self, topic):
            expect(topic).to_be_instance_of(Gamess)

        def should_have_a_tempdir(self, topic):
            #expect(topic.tempdir).to_be_str()
            expect(isinstance(topic.tempdir, str)).to_be_true()
