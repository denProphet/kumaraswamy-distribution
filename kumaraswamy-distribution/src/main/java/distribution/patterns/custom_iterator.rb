class CustomIterator

  #конструктор
  def initialize(array)
    @array = array
    @index = 0
  end

  #вернуть первый елемент
  def first
    @array[0]
  end

#возвращать следующий елемент
#увеличивать счетчик +1
  def next
    value = @array[@index]
    @index += 1
    value
  end

#возвращает true, если индекс итератора указывает на последний элемент в коллекции
  def is_done?
    @index == @array.length - 1
  end

#елемент по текущему индексу
  def current_item
    @array[@index]
  end

#определен ли следующий елемент?
  def has_next?
    @index < @array.length
  end

  #возвращать предыдущий елемент
  #уменьшать счетчик -1
  def previous
    value = @array[@index]
    @index -= 1
    value
  end

#определен ли предыдущий елемент?
  def has_previous?
    @index >= 0
  end

#сбросить счетчик
  def rewind
    @index = 0
  end
end
