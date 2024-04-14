import jg_impl_class as oc
import jg_impl_class_simplified as ocs

user_input = {}

print('--------- Normal ----------')

normal_cascade = oc.OxygenCascade(user_input)
print(normal_cascade.alt, normal_cascade.RR)
normal_cascade_out = normal_cascade.run_model()
print(normal_cascade_out)

user_input['alt'] = [8200, 'm']
user_cascade = oc.OxygenCascade(user_input)
print(user_cascade.alt, user_cascade.RR)
user_cascade_out = user_cascade.run_model()
print(user_cascade_out)

print('--------- Simplified ----------')

normal_cascade = ocs.OxygenCascadeSimplified({})
print(normal_cascade.alt, normal_cascade.RR)
normal_cascade_out = normal_cascade.run_simplified_model()
print(normal_cascade_out)

user_input['alt'] = [8848, 'm']
user_cascade = ocs.OxygenCascadeSimplified(user_input)
print(user_cascade.alt, user_cascade.RR)
user_cascade_out = user_cascade.run_simplified_model()
print(user_cascade_out)
