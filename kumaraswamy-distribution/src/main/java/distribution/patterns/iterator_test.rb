Pizzas_names = %w(Neapolitan Chicago Sicilian)


def external_iterator_impl_test
require_relative "custom_iterator.rb"

pizza_iterator = CustomIterator.new(Pizzas_names)

puts 'Кастомный итератор. По возврастанию '
while pizza_iterator.has_next?
  puts pizza_iterator.next
end

print 'Текущий индекс: '
puts(pizza_iterator.instance_variable_get( "@index" ))

puts 'Кастомный итератор. По убыванию '
while pizza_iterator.has_previous?
  puts pizza_iterator.previous
end

pizza_iterator.rewind

print 'Текущий индекс: '
puts(pizza_iterator.instance_variable_get( "@index" ))

puts 'Кастомный итератор. По возврастанию '
while pizza_iterator.has_next?
  puts pizza_iterator.next
end
print 'Текущий индекс: '
puts(pizza_iterator.instance_variable_get( "@index" ))

end

def internal_iterator_test
  puts 'Внутренний each итератор. По возврастанию '
  Pizzas_names.each do |name|
    puts name
  end
end

external_iterator_impl_test

internal_iterator_test
